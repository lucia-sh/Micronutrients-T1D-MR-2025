
#steiger filtering (remove SNPs that explain more variance in the outcome than in the exposure)
library(TwoSampleMR)

# dat is the MR harmonised file


dat$samplesize.exposure<- #exposuresamplesize
dat$samplesize.outcome<-#outcomesamplesize 
dat$pval.outcome<-dat$p_value

#fill in these numbers 
dat$ncase.outcome<-
dat$ncontrol.outcome<-
dat$prevalence.outcome<- ("n of cases"/("samplesize"))

dat$effective_n.exposure<- dat$samplesize.exposure
dat$effective_n.outcome <- 4/(1/"n of cases" + 1/"n of controls")


#calculate the variance explained 

dat$rsq.exposure <- (get_r_from_bsen(dat$beta.exposure,dat$se.exposure,241))^2 

dat$rsq.outcome <- get_r_from_lor(
  lor = (dat$beta.outcome),
  af = dat$eaf.outcome,
  ncase = dat$ncase.outcome,
  ncontrol = dat$ncontrol.outcome,
  prevalence = dat$prevalence.outcome,
  model = "logit",
  correction = FALSE)^2



# Apply Steiger filtering
mr_steig_filt <- steiger_filtering(dat)
mr_steig_filt2<-mr_steig_filt[,c('SNP','id.outcome','id.exposure','steiger_dir','steiger_pval')]