#-------------------*
# load packages
#-------------------*
source ("R/packages.R")

#-------------------*
# load data
#-------------------*

dados <- read.csv (here("data","data_anurans.csv"),h=T,sep=";")
dados <- dados [-674,]# the last line has NAs

#-------------------*
# formatting data
#-------------------*

dados$DATA <-as.Date (dados$DATA, "%d/%m/%Y")
# defining seasons
dados <- cbind (dados, ESTACAO=NA)
dados$ESTACAO <- ifelse (dados$DATA >= "2017-12-21" & dados$DATA <= "2018-03-20", "summer",dados$ESTACAO)
dados$ESTACAO <- ifelse (dados$DATA >= "2017-09-22" & dados$DATA <= "2017-12-21", "spring",dados$ESTACAO)
dados$ESTACAO <- ifelse (dados$DATA >= "2017-06-21" & dados$DATA <= "2017-09-22", "winter",dados$ESTACAO)
dados$ESTACAO <- ifelse (dados$DATA >= "2017-03-20" & dados$DATA <= "2017-06-21", "autumn",dados$ESTACAO)

# defining the month of sampling
dados <- cbind (dados, MONTH = format(dados$DATA,"%B"))

#-----------------------------------------------#
# bind a column with the basic sampling unit
# combination of season, month, area, habitat, transection 
#-----------------------------------------------#

dados <- cbind (dados, interaction=with (dados, interaction (ESTACAO,MONTH,AREA,HABITAT,TRANSECCAO)))

# composition matrix
matriz_sp_mes <- cast (interaction+ESTACAO+MONTH+AREA+HABITAT+TRANSECCAO+LATITUDE+LONGITUDE ~ ESPECIE,data=dados,
                       fun.aggregate=sum,value = "N_IND",na.rm=T)

# removing colname == NA
matriz_sp_mes<- matriz_sp_mes[,-which(colnames(matriz_sp_mes) == "NA")]

# editing the 'HABITAT' factor
matriz_sp_mes$HABITAT [which (matriz_sp_mes$HABITAT == "Nativa")] <- "Natural"
matriz_sp_mes$HABITAT <- factor(matriz_sp_mes$HABITAT,
                                levels = c("Natural",
                                           "Araucaria",
                                           "Pinus"))

# binding observed richness (counting of unique species)
frog_species <- unique (dados$ESPECIE)[-which(is.na(unique (dados$ESPECIE)))]

# numero observado de especies
matriz_sp_mes <- cbind  (matriz_sp_mes, OBS_RIQ=rowSums(decostand (matriz_sp_mes [,which(colnames(matriz_sp_mes) %in% frog_species)],
                                                                   "pa"))) 
# contagem de individuos
matriz_sp_mes <- cbind  (matriz_sp_mes, CONTAGEM=rowSums(matriz_sp_mes [,which(colnames(matriz_sp_mes) %in% frog_species)])) 
# plot( matriz_sp_mes$CONTAGEM,matriz_sp_mes$OBS_RIQ)

# indice de equidade na dist da abundancia
matriz_sp_mes<- cbind (matriz_sp_mes,
                       SIMPSON= diversity(matriz_sp_mes [,which(colnames(matriz_sp_mes) %in% frog_species)],
                                          index = "simpson", MARGIN = 1, base = exp(1)))

## A unica maneira de "livrar" o numero de especies do numero de individuos ?
## atraves de uma regressao entre riq ~ ind, e  entao usar os residuos

matriz_sp_mes <- cbind (matriz_sp_mes,
                        RES_RIQ= lm (OBS_RIQ ~ CONTAGEM, data=matriz_sp_mes)$res)

#---------------------------------#
# SPECIES COMPOSITION ANALYSIS
# NMDS
#---------------------------------#

composicao <-decostand (matriz_sp_mes[,which(colnames(matriz_sp_mes) %in% frog_species)],"hell")
composicao <- composicao[,which (colSums (matriz_sp_mes [,which(colnames(matriz_sp_mes) %in% frog_species)]) >=3)]

# NMDS with one axis to analyze composition diff
# the other with a single dimension (k = 1) to use the scores as a response variable in the GLMMs
species_composition_one_axis <- metaMDS (composicao,"euclidian",
                                         k=1,
                                         try = 20, 
                                         trymax = 20)
# goodness of fit
png (here("output","fit_NMDS.png"))
stressplot (species_composition_one_axis)
dev.off()

#  (k = 2) to visualize the clustering patterns formed by the distribution of sampling units
species_composition_two_axes <- metaMDS (composicao,"euclidean",
                                         k=2,
                                         try = 20, 
                                         trymax = 200)
# goodness of fit
stressplot (species_composition_two_axes)


# SAVE PLOT
pdf(file=here ("output","vectorized","species_composition_one_dim.pdf"),width=5, height=30,family="serif")

plot(species_composition_one_axis, 
     choices = 1, 
     display = c("sp","sites"),#c("sp", "wa", "cn"),
     ylim=c(-1,1),xlim=c(-1,1),
     cex=0)
dev.off()

# SAVE PLOT
pdf(file=here ("output","vectorized","species_composition_two_dim.pdf"),width=5, height=5,family="serif")

plot(species_composition_two_axes, 
     choices = c(1, 2), 
     display = "sites",#c("sp", "wa", "cn"),
     ylim=c(-1,1),xlim=c(-1,1),
     cex=0)

# point type to plot
tipos_pontos <- matriz_sp_mes$HABITAT
levels (tipos_pontos)[which(levels(tipos_pontos) == "Araucaria")] <- 1
levels (tipos_pontos)[which(levels(tipos_pontos) == "Pinus")] <- 0
levels (tipos_pontos)[which(levels(tipos_pontos) == "Natural")] <- 2
# plot points
points(species_composition_two_axes,
       pch=as.numeric(paste(tipos_pontos)),
       display = "sites",
       cex=1.2)
# text and legend
text (species_composition_two_axes,"species",cex=0.9)
legend("topleft", c("Natural","Araucaria","Pinus"), 
       col="black",pch=c(2,1,0),bty="n",cex=1)

text (2,3,"Stress=0.043",cex=1)

dev.off()

## COLANDO OS SCORES NA MATRIZ COM DADOS
matriz_sp_mes<- cbind(matriz_sp_mes, 
                      axis1=species_composition_one_axis$points)

# ------------------------ #
# SPATIAL AUTOCORRELATION
#------------------------------------------#
# help with autocorrelation
# he shows how to do by hand
# https://mgimond.github.io/Spatial/spatial-autocorrelation-in-r.html
#------------------------------------------#

# creating spatial points df
sp_points <- SpatialPoints(matriz_sp_mes[,c("LATITUDE","LONGITUDE")],
                           proj4string=CRS("+proj=longlat +datum=WGS84"))

# creating spatial neighborhood
nb <- knn2nb(knearneigh(sp_points))
# weights
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
neighwb <- nb2WB(nb)

# test
(MT.rich <- moran.mc(matriz_sp_mes$OBS_RIQ,lw,nsim=999))
(MT.cont<-moran.mc(matriz_sp_mes$CONTAGEM,lw,nsim=999))
(MT.comp<-moran.mc(matriz_sp_mes$MDS1,lw,nsim=999))
(MT.cuvieri<-moran.mc(matriz_sp_mes$Physalaemus_cuvieri,lw,nsim=999))
(MT.gracilis<-moran.mc(matriz_sp_mes$Physalaemus_gracilis,lw,nsim=999))
(MT.henseli<-moran.mc(matriz_sp_mes$Rhinella_henseli,lw,nsim=999))
(MT.araucaria<-moran.mc(matriz_sp_mes$Adenomera_araucaria,lw,nsim=999))

plot(MT.rich, main="", las=1)
plot(MT.cont, main="", las=1)
plot(MT.comp, main="", las=1)

# yes, we have spatial autocorrelation
# defining neighborhood
require(rgeos)
sp_points_buffer <- gBuffer(sp_points,byid=T,width=1)

# neighborhood between polygons
neigh <- poly2nb(sp_points_buffer)
neighwb <- nb2WB(neigh)

## spatial structure
num <- neighwb$num
#num <- rep (nvizinhos,496)## n neighbors
#adj<- unlist(neigh$neighbours)
W <- matrix (NA, nrow =length (num), 
             ncol=length (num)) ## mmatrix of sites (rows) by neighbors (cols)
# adjacency matrix
adj_matrix <- matrix(neighwb$adj,nrow = dim(W), ncol=25,byrow=T)
W_adj <- do.call (rbind, 
                  
                  lapply (seq (1,nrow(W)), function (i) {
  
                  W [i,adj_matrix [i,]]<-1 ;
                  W [i,]

  
}))

W_adj [is.na(W_adj)] <- 0
# ------------------- #
# mixed models
# with spatial autocorrelation
# and area nested in season
# help 
# http://www.flutterbys.com.au/stats/tut/tut9.1.html
# https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html
# ------------------- #

# MCMC settings
ni <- 50000
nb <- 40000
nt <- 20

# counting of individuals
# RANEF: area nested into season
m1_cont <- brm(CONTAGEM~HABITAT+
                 (1|AREA/ESTACAO) +
                 car(W),
          data=matriz_sp_mes,
          data2 = list(W = W_adj),
          family=poisson(link="log"),
          chains = 3,
          iter = ni,
          warmup = nb,
          thin = nt)

# counting of individuals
# RANEF: season
# GROUPING: AREA
m2_cont <- brm(CONTAGEM~HABITAT+
                 (1|ESTACAO)+
                 car(W),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m3_cont <- brm(CONTAGEM~HABITAT+
                 (1|AREA)+
                 car(W),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

# alternative models
m1_cont_noCAR <- brm((CONTAGEM)~HABITAT+
                 (1|AREA/ESTACAO),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m2_cont_noCAR <- brm(CONTAGEM~HABITAT+
                       (1|ESTACAO),
                     data=matriz_sp_mes,
                     data2 = list(W = W_adj),
                     family=poisson(link="log"),
                     chains = 3,
                     iter = ni,
                     warmup = nb,
                     thin = nt)
pairs(m2_cont_noCAR)

m3_cont_noCAR <- brm(CONTAGEM~HABITAT+
                       (1|AREA),
                     data=matriz_sp_mes,
                     data2 = list(W = W_adj),
                     family=poisson(link="log"),
                     chains = 3,
                     iter = ni,
                     warmup = nb,
                     thin = nt)


save (m1_cont,m2_cont,m3_cont, file=here("output","brms_counting.RData"))
save (m1_cont_noCAR,m2_cont_noCAR,m3_cont_noCAR, file=here("output","brms_counting_m4.RData"))

# richness
m1_rich <- brm(OBS_RIQ~HABITAT+
                 (1|AREA/ESTACAO) +
                 car(W),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)


# richness
m2_rich <- brm(OBS_RIQ~HABITAT+
                 (1|ESTACAO)+
                 car(W),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m3_rich <- brm(OBS_RIQ~HABITAT+
                 (1|AREA)+
                 car(W),
               data=matriz_sp_mes,
               data2 = list(W = W_adj),
               family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

save (m1_rich,m2_rich,m3_rich, file=here("output","brms_richness.RData"))

# speciesi composition
#
m1_comp <- brm(MDS1 ~ HABITAT+
                 (1|AREA/ESTACAO),
               data=matriz_sp_mes,
               family=gaussian(link="identity"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m2_comp <- brm(MDS1 ~ HABITAT+
                 (1|ESTACAO),
               data=matriz_sp_mes,
               family=gaussian(link="identity"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)


m3_comp <- brm(MDS1 ~ HABITAT+
                 (1|AREA),
               data=matriz_sp_mes,
               family=gaussian(link="identity"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

save (m1_comp,m2_comp,m3_comp, file=here("output","brms_composition.RData"))

# cuvieri

m1_cuvieri <- brm(Physalaemus_cuvieri ~ HABITAT+
                 (1|AREA/ESTACAO)+
                 car(W),
                 data=matriz_sp_mes,
                 data2 = list(W = W_adj),
                 family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m2_cuvieri <- brm(Physalaemus_cuvieri ~ HABITAT+
                 (1|ESTACAO)+
                 car(W),
                 data=matriz_sp_mes,
                 data2 = list(W = W_adj),
                 family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

m3_cuvieri <- brm(Physalaemus_cuvieri ~ HABITAT+
                 (1|AREA)+
                   car(W),
                 data=matriz_sp_mes,
                 data2 = list(W = W_adj),
                 family=poisson(link="log"),
               chains = 3,
               iter = ni,
               warmup = nb,
               thin = nt)

save (m1_cuvieri,m2_cuvieri,m3_cuvieri, file=here("output","brms_cuvieri.RData"))


# alternative models
m1_cuvieri_noCAR <- brm(Physalaemus_cuvieri~HABITAT+
                       (1|AREA/ESTACAO),
                     data=matriz_sp_mes,
                     data2 = list(W = W_adj),
                     family=poisson(link="log"),
                     chains = 3,
                     iter = ni,
                     warmup = nb,
                     thin = nt)

m2_cuvieri_noCAR <- brm(Physalaemus_cuvieri~HABITAT+
                       (1|ESTACAO),
                     data=matriz_sp_mes,
                     data2 = list(W = W_adj),
                     family=poisson(link="log"),
                     chains = 3,
                     iter = ni,
                     warmup = nb,
                     thin = nt)

m3_cuvieri_noCAR <- brm(Physalaemus_cuvieri~HABITAT+
                       (1|AREA),
                     data=matriz_sp_mes,
                     data2 = list(W = W_adj),
                     family=poisson(link="log"),
                     chains = 3,
                     iter = ni,
                     warmup = nb,
                     thin = nt)


save (m1_cuvieri_noCAR,m2_cuvieri_noCAR,m3_cuvieri_noCAR, file=here("output","brms_cuvieri_alternative.RData"))

# gracilis 

m1_gracilis <- brm(Physalaemus_gracilis ~ HABITAT+
                    (1|AREA/ESTACAO)+
                     car(W),
                   data=matriz_sp_mes,
                   data2 = list(W = W_adj),
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

m2_gracilis <- brm(Physalaemus_gracilis ~ HABITAT+
                    (1|ESTACAO)+
                     car(W),
                   data=matriz_sp_mes,
                   data2 = list(W = W_adj),
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)


m3_gracilis <- brm(Physalaemus_gracilis ~ HABITAT+
                    (1|AREA)+
                     car(W),
                   data=matriz_sp_mes,
                   data2 = list(W = W_adj),
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

save (m1_gracilis,m2_gracilis,m3_gracilis, file=here("output","brms_gracilis.RData"))

# alternative


# alternative models
m1_gracilis_noCAR <- brm(Physalaemus_gracilis~HABITAT+
                          (1|AREA/ESTACAO),
                        data=matriz_sp_mes,
                        data2 = list(W = W_adj),
                        family=poisson(link="log"),
                        chains = 3,
                        iter = ni,
                        warmup = nb,
                        thin = nt,
                        moment_match = TRUE)

m2_gracilis_noCAR <- brm(Physalaemus_gracilis~HABITAT+
                          (1|ESTACAO),
                        data=matriz_sp_mes,
                        data2 = list(W = W_adj),
                        family=poisson(link="log"),
                        chains = 3,
                        iter = ni,
                        warmup = nb,
                        thin = nt,
                        moment_match = TRUE)

m3_gracilis_noCAR <- brm(Physalaemus_gracilis~HABITAT+
                          (1|AREA),
                        data=matriz_sp_mes,
                        data2 = list(W = W_adj),
                        family=poisson(link="log"),
                        chains = 3,
                        iter = ni,
                        warmup = nb,
                        thin = nt,
                        moment_match = TRUE)


save (m1_gracilis_noCAR,m2_gracilis_noCAR,m3_gracilis_noCAR, file=here("output","brms_gracilis_alternative.RData"))


# R henseli

m1_henseli <- brm(Rhinella_henseli ~ HABITAT+
                     (1|AREA/ESTACAO),
                   data=matriz_sp_mes,
                   family=poisson(link="log"),
                   chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

m2_henseli <- brm(Rhinella_henseli ~ HABITAT+
                     (1|ESTACAO),
                   data=matriz_sp_mes,
                   family=poisson(link="log"),
                   chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)


m3_henseli <- brm(Rhinella_henseli ~ HABITAT+
                     (1|AREA),
                   data=matriz_sp_mes,
                   family=poisson(link="log"),
                   chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

save (m1_henseli,m2_henseli,m3_henseli, file=here("output","brms_henseli.RData"))

# adenomera araucaria

m1_araucaria <- brm(Adenomera_araucaria ~ HABITAT+
                    (1|AREA/ESTACAO),
                  data=matriz_sp_mes,
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

m2_araucaria <- brm(Adenomera_araucaria ~ HABITAT+
                    (1|ESTACAO),
                  data=matriz_sp_mes,
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

m3_araucaria <- brm(Adenomera_araucaria ~ HABITAT+
                    (1|AREA),
                  data=matriz_sp_mes,
                  family=poisson(link="log"),
                  chains = 3,
                  iter = ni,
                  warmup = nb,
                  thin = nt)

save (m1_araucaria,m2_araucaria,m3_araucaria, file=here("output","brms_araucaria.RData"))
