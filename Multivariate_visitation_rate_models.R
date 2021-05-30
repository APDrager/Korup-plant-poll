
#R script for running analysis used to estimate visitation rate parameters in the three models in Drager et al. 
#'What structures diurnal visitation rates to flowering trees in an Afrotropical lowland rain-forest understory?' 

#Author: A. P. Drager

#Helpful reference: #https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html
#Thanks also to this reference for much coding help: https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/

library(dplyr)
library(ggplot2)
library(brms)
library(bayesplot)
library(tidybayes)
library(tidyr)
library(broom)
library(ggcorrplot)
library(phytools)
library(MCMCglmm)


dat<-read.csv("Korup-floral-visitors_rawdata.csv") #input plot dataframe 'dat'

#Model Priors: biologically relevant prior avoids extreme (impossible) values for standard deviations of visitation rates. 
#all other parameters use brms default priors
eprior <-set_prior("exponential(1)", class = "sd", resp = c("ants","bees","beetles","bugs","neodipt","dipt"))

## Multivariate models:

#Model 1. Full model--scent categories & day

bfd_ants <- bf(ants ~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bfd_bees <- bf(bees ~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bfd_beetles <- bf(beetles ~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bfd_bugs <- bf(bugs ~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bfd_neo_dipt <- bf(neo_dipt~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bfd_dipt <- bf(dipt~ 0 + fermented + unscented + sweet + day + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))


Model1 <- brm(bfd_ants + bfd_bees + bfd_beetles + bfd_bugs + bfd_neo_dipt + bfd_dipt,
                      data = dat,  prior=eprior, family = zero_inflated_poisson ,  control = list(adapt_delta = 0.99),
                      chains=2, cores=2, iter = 10000, warmup = 3000)

#Save model fit to avoid running again
save(Model1, file = "Model1.rda")

#View detailed posterior mean estimates
posterior_summary(Model1) %>%round(digits = 2)

#Model 2. Scent category model

bf_ants <- bf(ants ~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bf_bees <- bf(bees ~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bf_beetles <- bf(beetles ~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bf_bugs <- bf(bugs ~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bf_neo_dipt <- bf(neo_dipt~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))
bf_dipt <- bf(dipt~ 0 + fermented + unscented + sweet + offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code))

Model2 <- brm(bf_ants + bf_bees + bf_beetles + bf_bugs + bf_neo_dipt + bf_dipt,
                   data = dat,  prior=eprior, family = zero_inflated_poisson ,  control = list(adapt_delta = 0.99),
                   chains=4, iter = 10000, warmup = 3000)

save(Model2, file = "Model2.rda")

#Model 2 with phylogenetically-structured species random effects:

phylo<-read.nexus("Korup_PP_phylo.nex") 

inv.phylo <- MCMCglmm::inverseA(phylo, nodes = "TIPS", scale = TRUE) 
A <- solve(inv.phylo$Ainv)
colnames(A) <- rownames(A) <- rownames(inv.phylo$Ainv) #matrix with all the tips, not just the ones in the dataset...huh, why was it like this for Michael vs just having our 23 species?


Model2_phylo <- brm(bf_ants + bf_bees + bf_beetles + bf_bugs + bf_neo_dipt + bf_dipt,
                    data = dat,  prior=eprior, family = zero_inflated_poisson ,  control = list(adapt_delta = 0.99),
                    cov_ranef=list(code=A),
                    chains=4, iter = 10000, warmup = 3000)

save(Model2_phylo, file = "Model2_phylo.rda")

#Model 3. Random effects only model:

Model3<- brm( mvbind(ants, bees, beetles, bugs, neo_dipt, dipt) ~ 1 +offset(log(ob_time_min/60)) + (1|q|tree_tag) + (1|p|code),
                 data = dat,  prior=eprior, family = zero_inflated_poisson ,  control = list(adapt_delta = 0.99),
                 chains=4, iter = 10000, warmup = 3000)

save(Model3, file = "Model3.rda")


#Examine model performance using deviance metric (ratio of N effective samples to N samples):
neff_ratio(Model1) %>% 
  mcmc_neff()

#Comparing model fits va 10-fold random cross validation:

kfold(Model1, Model2, Model2_phylo, Model3)


#Save posterior samples for future use:
postModel1 <- posterior_samples(Model1, add_chain = T)
postModel2 <- posterior_samples(Model2, add_chain = T)
postModel3 <- posterior_samples(Model3, add_chain = T)


#Pairwise tests for significant differences in zero-inflation rates between groups
postModel1 %>%
  transmute(dif_abz = zi_ants - zi_bees,
            dif_acz = zi_ants - zi_beetles,
            dif_ahz = zi_ants - zi_bugs,
            dif_anz = zi_ants - zi_neodipt,
            dif_adz = zi_ants - zi_dipt,
            dif_bcz = zi_bees - zi_beetles,
            dif_bhz = zi_bees - zi_bugs,
            dif_bnz = zi_bees - zi_neodipt,
            dif_bdz = zi_bees - zi_dipt,
            dif_chz = zi_beetles - zi_bugs,
            dif_cnz = zi_beetles - zi_neodipt,
            dif_cdz = zi_beetles- zi_dipt,
            dif_hnz = zi_bugs - zi_neodipt,
            dif_hdz = zi_bugs - zi_dipt,
            dif_ndz = zi_neodipt - zi_dipt)%>%
  gather() %>%
  mutate(key = factor(key, levels = c("dif_abz", "dif_acz", "dif_ahz","dif_anz", "dif_adz", "dif_bcz","dif_bhz", "dif_bnz", "dif_bdz", "dif_chz", "dif_cnz", "dif_cdz", "dif_hnz", "dif_hdz", "dif_ndz"))) %>%
  group_by(key) %>%
  median_qi(value, .width = .95) %>% 
  mutate_if(is.double, round, digits = 2)-> VisRate_ZI_Pairwise_tests_median_qi90e

#Pairwise tests for significant difference in visitation rates WITHIN visitor groups to different scent categories
postModel2 %>%
  transmute(dif_afu = b_ants_fermented - b_ants_unscented,
            dif_afs = b_ants_fermented - b_ants_sweet,
            dif_aus = b_ants_unscented - b_ants_sweet,
            dif_bfu = b_bees_fermented - b_bees_unscented,
            dif_bfs = b_bees_fermented - b_bees_sweet,
            dif_bus = b_bees_unscented - b_bees_sweet,
            dif_cfu = b_beetles_fermented - b_beetles_unscented,
            dif_cfs = b_beetles_fermented - b_beetles_sweet,
            dif_cus = b_beetles_unscented - b_beetles_sweet,
            dif_hfu = b_bugs_fermented - b_bugs_unscented,
            dif_hfs = b_bugs_fermented - b_bugs_sweet,
            dif_hus = b_bugs_unscented - b_bugs_sweet,
            dif_nfu = b_neodipt_fermented - b_neodipt_unscented,
            dif_nfs = b_neodipt_fermented - b_neodipt_sweet,
            dif_nus = b_neodipt_unscented - b_neodipt_sweet,
            dif_dfu = b_dipt_fermented - b_dipt_unscented,
            dif_dfs = b_dipt_fermented - b_dipt_sweet,
            dif_dus = b_dipt_unscented - b_dipt_sweet)%>%
  gather() %>%
  mutate(key = factor(key, levels = c("dif_afu", "dif_afs", "dif_aus","dif_bfu", "dif_bfs", "dif_bus","dif_cfu", "dif_cfs", "dif_cus", "dif_hfu", "dif_hfs", "dif_hus", "dif_nfu", "dif_nfs", "dif_nus" ,"dif_dfu", "dif_dfs", "dif_dus"))) %>%
  group_by(key) %>%
  median_qi(value, .width = .90) %>% 
  mutate_if(is.double, round, digits = 2)-> VisRate_Scent_WithinVisGrp_pairtests

#Pairwise tests for significant difference in visitation rates AMONG visitor groups to different scent categories
postModel2 %>%
  transmute(dif_abf = b_ants_fermented - b_bees_fermented,
            dif_abu = b_ants_unscented - b_bees_unscented,
            dif_abs = b_ants_sweet - b_bees_sweet,
            dif_acf = b_ants_fermented - b_beetles_fermented,
            dif_acu = b_ants_unscented - b_beetles_unscented,
            dif_acs = b_ants_sweet - b_beetles_sweet,
            dif_ahf = b_ants_fermented - b_bugs_fermented,
            dif_ahu = b_ants_unscented - b_bugs_unscented,
            dif_ahs = b_ants_sweet - b_bugs_sweet,
            dif_anf = b_ants_fermented - b_neodipt_fermented,
            dif_anu = b_ants_unscented - b_neodipt_unscented,
            dif_ans = b_ants_sweet - b_neodipt_sweet,
            dif_adf = b_ants_fermented - b_dipt_fermented,
            dif_adu = b_ants_unscented - b_dipt_unscented,
            dif_ads = b_ants_sweet - b_dipt_sweet,
            dif_bcf = b_bees_fermented - b_beetles_fermented,
            dif_bcu = b_bees_unscented - b_beetles_unscented,
            dif_bcs = b_bees_sweet - b_beetles_sweet,
            dif_bhf = b_bees_fermented - b_bugs_fermented,
            dif_bhu = b_bees_unscented - b_bugs_unscented,
            dif_bhs = b_bees_sweet - b_bugs_sweet,
            dif_bnf = b_bees_fermented - b_neodipt_fermented,
            dif_bnu = b_bees_unscented - b_neodipt_unscented,
            dif_bns = b_bees_sweet - b_neodipt_sweet,
            dif_bdf = b_bees_fermented - b_dipt_fermented,
            dif_bdu = b_bees_unscented - b_dipt_unscented,
            dif_bds = b_bees_sweet - b_dipt_sweet,
            dif_chf = b_beetles_fermented - b_bugs_fermented,
            dif_chu = b_beetles_unscented - b_bugs_unscented,
            dif_chs = b_beetles_sweet - b_bugs_sweet,
            dif_cnf = b_beetles_fermented - b_neodipt_fermented,
            dif_cnu = b_beetles_unscented - b_neodipt_unscented,
            dif_cns = b_beetles_sweet - b_neodipt_sweet,
            dif_cdf = b_beetles_fermented - b_dipt_fermented,
            dif_cdu = b_beetles_unscented - b_dipt_unscented,
            dif_cds = b_beetles_sweet - b_dipt_sweet,
            dif_hnf = b_bugs_fermented - b_neodipt_fermented,
            dif_hnu = b_bugs_unscented - b_neodipt_unscented,
            dif_hns = b_bugs_sweet - b_neodipt_sweet,
            dif_hdf = b_bugs_fermented - b_dipt_fermented,
            dif_hdu = b_bugs_unscented - b_dipt_unscented,
            dif_hds = b_bugs_sweet - b_dipt_sweet,
            dif_ndf = b_neodipt_fermented - b_dipt_fermented,
            dif_ndu = b_neodipt_unscented - b_dipt_unscented,
            dif_nds = b_neodipt_sweet - b_dipt_sweet) %>%
 gather() %>%
 mutate(key = factor(key, levels = c("dif_abf", "dif_abu", "dif_abs","dif_acf", "dif_acu", "dif_acs","dif_ahf", "dif_ahu", "dif_ahs", "dif_anf", "dif_anu", "dif_ans", "dif_adf", "dif_adu", "dif_ads" ,"dif_bcf", "dif_bcu", "dif_bcs","dif_bhf", "dif_bhu", "dif_bhs", 
                                      "dif_bnf", "dif_bnu", "dif_bns","dif_bdf", "dif_bdu", "dif_bds", "dif_chf", "dif_chu", "dif_chs","dif_cnf", "dif_cnu", "dif_cns", "dif_cdf", "dif_cdu", "dif_cds", "dif_hnf", "dif_hnu", "dif_hns", "dif_hdf", "dif_hdu", "dif_hds","dif_ndf", "dif_ndu", "dif_nds"))) %>%
 group_by(key) %>%
 median_qi(value, .width = .90) %>% 
 mutate_if(is.double, round, digits = 2)-> VisRate_Scent_AmongVisGrps_pairtests


#Example of how we can obtain rates from parameter estimates: 
postModel3 %>%
  select(starts_with("b_")) %>% 
  gather() %>% # this is how we might add the grand mean to the tree-species deviations
  mutate(value = exp(value)) %>% 
  group_by(key) %>%
  median_qi(value, .width = c(.95)) %>% 
  mutate_if(is.double, round, digits = 2) -> beta_rates


#Making the plots for correlation in insect group visitation rates to tree spp. and to individual trees (tree tags)
#Manuscript Figure S3.
#make matrix for correlation in tree spp: 
tidy(Model3)->Model3tidy 
Model3tidy %>% slice(19:33)->cor_code  #just correlation estimates for rates to tree spp.
cor_code$estimate->corrcode  
corcode<- matrix(0, nrow = 6, ncol = 6) #create matrix from vector:
mindex <- matrix(1:36, nrow = 6, ncol = 6)
corcode[mindex[upper.tri(mindex)]] <- corrcode
corcode[lower.tri(corcode)] <- t(corcode)[lower.tri(corcode)]
diag(corcode) <- 1
corcode
colnames(corcode) <- c("Ants", "Bees", "Beetles", "Bugs", "Delicate flies", "Robust flies")
rownames(corcode) <- c("Ants", "Bees", "Beetles", "Bugs", "Delicate flies", "Robust flies")
corr<-round(corcode, 2)
#plot matrix
ggcorrplot(corr,method="circle", type="upper", title="Tree species  \n", lab=T, lab_size = 3)->code_corrplot


#make matrix for correlation in tree tags: 
Model3tidy %>% slice(34:48)->cor_tag  
cor_tag$estimate->corrtag  
cortag<- matrix(0, nrow = 6, ncol = 6)
mindex <- matrix(1:36, nrow = 6, ncol = 6)
cortag[mindex[upper.tri(mindex)]] <- corrtag
cortag[lower.tri(cortag)] <- t(cortag)[lower.tri(cortag)]
diag(cortag) <- 1
cortag
colnames(cortag) <- c("Ants", "Bees", "Beetles", "Bugs", "Delicate flies", "Robust flies")
rownames(cortag) <- c("Ants", "Bees", "Beetles", "Bugs", "Delicate flies", "Robust flies")
cort<-round(cortag, 2)
ggcorrplot(cort,method="circle", type="upper", title= "Individual trees \n", lab=T, lab_size = 3)-> tag_corrplot

grid.arrange(code_corrplot, tag_corrplot, ncol=2)

