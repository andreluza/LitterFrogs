## ----------------------------------------------
# CODE FOR INTERPRETATION

#-------------------*
# load packages
#-------------------*
source ("R/packages.R")

#-------------------*
# load results
#-------------------*

load (here("output","brms_richness.RData"))
load (here("output","brms_counting.RData"))
load (here("output","brms_counting_m4.RData"))
load (here("output","brms_composition.RData"))
load (here("output","brms_cuvieri.RData"))
load (here("output","brms_cuvieri_alternative.RData"))
load (here("output","brms_gracilis.RData"))
load (here("output","brms_gracilis_alternative.RData"))
load (here("output","brms_henseli.RData"))
load (here("output","brms_araucaria.RData"))

# check model fit
# richness 
loo(m1_rich)
loo(m2_rich)
loo(m3_rich)

# abundance 
loo(m1_cont)
loo(m2_cont)
loo(m3_cont)

# alternative of abundance
loo(m1_cont_noCAR)
loo(m2_cont_noCAR)
loo(m3_cont_noCAR)


# composition 
loo(m1_comp)
loo(m2_comp)
loo(m3_comp)

# cuvieri
loo(m1_cuvieri)
loo(m2_cuvieri)
loo(m3_cuvieri)

# cuvieri
loo(m1_cuvieri_noCAR)
loo(m2_cuvieri_noCAR)
loo(m3_cuvieri_noCAR)


# gracilis
loo(m1_gracilis)
loo(m2_gracilis)
loo(m3_gracilis)


# alternative
loo(m1_gracilis_noCAR)
loo(m2_gracilis_noCAR)
loo(m3_gracilis_noCAR)

# henseli
loo(m1_henseli)
loo(m2_henseli)
loo(m3_henseli)

# araucaria
loo(m1_araucaria)
loo(m2_araucaria)
loo(m3_araucaria)

#pp_check(m1_rich,nsamples=100)
#pp_check(m2_rich,nsamples=100)
#pp_check(m3_rich,nsamples=100)

# generate a summary of the results
summary(m2_rich)

# plot the MCMC chains as well as the posterior distributions
plot(m2_rich, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_rich))

# plot conditional effects for each predictor
plot((m2_rich))

pdf (here("output","vectorized", "richness.pdf"),width=3,height=3)

plot(conditional_effects(m2_rich,
                    method="fitted",
                    re_formula=NA,
                    robust=F),
                    
     theme = theme_classic(),
     points=F) 

dev.off()

## total abundance

# generate a summary of the results
summary(m2_cont_noCAR)

# plot the MCMC chains as well as the posterior distributions
plot(m2_cont_noCAR, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_cont_noCAR))

# plot conditional effects for each predictor
plot((m2_cont_noCAR))

pdf (here("output","vectorized", "abundance.pdf"),width=3,height=3)

plot(conditional_effects(m2_cont_noCAR,
                         method="fitted",
                         re_formula=NA,
                         robust=F),
     
     theme = theme_classic(),
     points=F) 

dev.off()

# composition

# generate a summary of the results
summary(m2_comp)

# plot the MCMC chains as well as the posterior distributions
plot(m2_comp, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_comp))

# plot conditional effects for each predictor
plot((m2_comp))

pdf (here("output","vectorized", "composition.pdf"),width=3,height=3)

plot(conditional_effects(m2_comp,
                         method="fitted",
                         re_formula=NA,
                         robust=F),
     
     theme = theme_classic(),
     points=F) 

dev.off()

# cuvieri

# generate a summary of the results
summary(m2_cuvieri_noCAR)

# plot the MCMC chains as well as the posterior distributions
plot(m2_cuvieri_noCAR, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_cuvieri_noCAR))

# plot conditional effects for each predictor
plot((m2_cuvieri_noCAR))

pdf (here("output","vectorized", "cuvieri.pdf"),width=3,height=3)

plot(conditional_effects(m2_cuvieri_noCAR,
                         method="fitted",
                         re_formula=NA,
                         robust=F),
     
     theme = theme_classic(),
     points=F) 

dev.off()


# carrizorum

# generate a summary of the results
summary(m2_gracilis_noCAR)

# plot the MCMC chains as well as the posterior distributions
plot(m2_gracilis_noCAR, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_gracilis_noCAR))

# plot conditional effects for each predictor
plot((m2_gracilis_noCAR))

pdf (here("output","vectorized", "carrizorum.pdf"),width=3,height=3)

plot(conditional_effects(m2_gracilis_noCAR,
                         method="fitted",
                         re_formula=NA,
                         robust=F),
     
     theme = theme_classic(),
     points=F) 

dev.off()


# henseli

# generate a summary of the results
summary(m2_henseli)

# plot the MCMC chains as well as the posterior distributions
plot(m2_henseli, ask = FALSE)

# predict responses based on the fitted model
head(predict(m2_henseli))

# plot conditional effects for each predictor
plot((m2_henseli))

pdf (here("output","vectorized", "henseli.pdf"),width=3,height=3)

plot(conditional_effects(m2_henseli,
                         method="fitted",
                         re_formula=NA,
                         robust=F),
     
     theme = theme_classic(),
     points=F) 

dev.off()
