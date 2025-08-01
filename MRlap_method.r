remotes::install_github("n-mounier/MRlap")

library(MRlap)
library(data.table)
library(dplyr)



# Define exposure and outcome names
exposure_name <- ""
outcome_name  <- "T1D"

exposure<-as.data.frame(fread('pathway to exposure'))
outcome<-as.data.frame(fread('pathway to outcome'))


#rename the columns of your exposure and outcome  dataframe as: CHR, POS, beta , se,  alt, ref, rsid, N

names(exposure)[names(exposure)=='chromosome']<-"CHR"
names(exposure)[names(exposure)=='base_pair_location']<-"POS"
names(exposure)[names(exposure)=='standard_error']<-"se"
names(exposure)[names(exposure)=='n']<-"N"
names(exposure)[names(exposure)=='effect_allele']<-"alt"
names(exposure)[names(exposure)=='other_allele']<-"ref"

names(outcome)[names(outcome)=="chromosome"]<-"CHR"
names(outcome)[names(outcome)=="base_pair_location"]<-"POS"
names(outcome)[names(outcome)=='standard_error'] <- 'N'
names(outcome)[names(outcome)=='Allele1'] <- 'alt'
names(outcome)[names(outcome)=='Allele2'] <- 'ref'
names(outcome)[names(outcome)=='SE'] <- 'se'
names(outcome)[names(outcome)=='Beta'] <- 'beta'


dat<-readRDS('the harmonized exposure file')

dat<-dat[dat$mr_keep==TRUE,]
exposure_filtered <- exposure %>% filter(rsid %in% dat$SNP)

result <- MRlap(
  exposure      = exposure,
  exposure_name = exposure_name,
  outcome       = outcome,
  outcome_name  = outcome_name,

# Set paths to required reference files for LD Score Regression (LDSC)
# 'ld' should point to the directory containing LD score files (e.g., eur_w_ld_chr)
# 'hm3' should point to the HapMap3 SNP list used to filter SNPs (e.g., w_hm3.snplist)

  ld            = "path to ld files eur_w_ld_chr",
  hm3           = "path to HapMap3 w_hm3.snplist",
  do_pruning    = FALSE,
  user_SNPsToKeep = exposure_filtered$rsid,
  MR_reverse    = 1e-3,
  verbose       = TRUE
)

#save MRLap file
saveRDS(result,'path')

#process_results
unlist(result[["MRcorrection"]])

# 'A' is the MRlap result object
MR_res <- result$MRcorrection

# Create a transposed one-row table with named columns
MR_table <- data.frame(
  Outcome              = 'T1D',
  Exposure             = '',
  Observed_Effect      = MR_res$observed_effect,
  Observed_SE          = MR_res$observed_effect_se,
  Observed_P           = MR_res$observed_effect_p,
  Corrected_Effect     = MR_res$corrected_effect,
  Corrected_SE         = MR_res$corrected_effect_se,
  Corrected_P          = MR_res$corrected_effect_p,
  Number_of_IVs        = MR_res$m_IVs,
  Test_Diff_Obs_vs_Corr= MR_res$test_difference,
  P_Diff_Obs_vs_Corr   = MR_res$p_difference,
  stringsAsFactors     = FALSE
)


print(MR_table)

LDSC_res <- result$LDSC

LDSC_table<-data.frame(
  Outcome              = 'T1D',
  Exposure             = '',
  h2_exp=LDSC_res$h2_exp,
  h2_exp_se=LDSC_res$h2_exp_se,
  int_exp=LDSC_res$int_exp,
  h2_out=LDSC_res$h2_out,
  h2_out_se=LDSC_res$h2_out_se,
  int_out=LDSC_res$int_out,
  gcov=LDSC_res$gcov,
  gcov_se=LDSC_res$gcov_se,
  rg=LDSC_res$rg,
  int_crosstrait=LDSC_res$int_crosstrait,
  int_cross_trait_se=LDSC_res$int_crosstrait_se,
  stringsAsFactors     = FALSE
)

print(LDSC_table)

GA_res <- result$GeneticArchitecture

GA_table <- data.frame(
  polygenicity=GA_res$polygenicity,
  perSNP_heritability=GA_res$perSNP_heritability
)


print(GA_table)
