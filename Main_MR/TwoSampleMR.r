
#R 

#Install packages 

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("MRCIEU/genetics.binaRies")


#Load Libraries

library(TwoSampleMR)
library(dplyr)
library(data.table)
library(genetics.binaRies)

#Upload our GWAS sumstats

freadf <- function(df){return(as.data.frame(fread(df)))}

exposure<-freadf('exposure pathway')
outcome<-freadf('outcome pathway')

#rename columns if needed
#Exposure Filtering

###P_value filtering 

exposure$p_value<-as.numeric(exposure$p_value)
exposure_test<-exposure[exposure$p_value<5e-6,]

exposure<-exposure_test 


###Remove rare variants 

exposure_test<-subset(exposure,effect_allele_frequency>0.01 & effect_allele_frequency<0.99)
exposure<-exposure_test 

###Remove HLA region

exposure_test<- exposure %>% filter(!(chromosome == 6 & base_pair_location >= 28510120 & base_pair_location <= 33480577))

exposure<-exposure_test 

###if we don't have frequency information

#we add 0.5 eaf if it is not available 
if (!"eaf.exposure" %in% names(exposure)) {
  exposure$eaf.exposure <- rep(0.5, nrow(exposure))
}

exposure<-exposure_test

###Remove Indels 

exposure_test<-exposure[nchar(exposure$effect_allele)== 1 &
                     nchar(exposure$other_allele)== 1 &
                        exposure$effect_allele %in% c("A" , "T", "C" , "G") &
                exposure$other_allele %in% c("A" ,"T", "C","G"),]

exposure<-exposure_test

outcome_test <- outcome[nchar(outcome$effect_allele) == 1 & 
                           nchar(outcome$other_allele) == 1 & 

                           outcome$effect_allele %in% c("A", "T", "C", "G") & 
                           outcome$other_allele %in% c("A", "T", "C", "G"), ]
outcome<-outcome_test

###Fstat Filtering 
exposure$N<-694649 #If no n (samplesize) col, set this value. Otherwise use the given n (sample_size)
exposure$R2<-get_r_from_bsen(exposure$beta, exposure$standard_error, exposure$n)^2
exposure$Fstat <- (exposure$R2/1)/((1-exposure$R2)/(exposure$n-1-1))

#Filtering 
exposure_test <- exposure[exposure$Fstat > 10,]

exposure<-exposure_test


###LD Clumping

#LD Clumping
#change columns as needed 

names(exposure)[names(exposure)=='chromosome']<-"chr_name"
names(exposure)[names(exposure)=='base_pair_location']<-"chrom_start"
names(exposure)[names(exposure)=='p_value']<-"pval.exposure"
names(exposure)[names(exposure)=='rsid']<-"SNP"

# LD Clumping
clumped_data <- data.frame()
for (i in 1:22) {
  cat("Clumping chromosome", i, "\n")
  df_chr <- exposure[exposure$chr_name == i, ]
  bfile_path <- paste0('/1000Genomes/wgs/plink/EUR/Chr', i, '.EUR') #pathway to 1000 Genomes EUR reference panel
  cat(bfile_path, "\n")
  
  # Wrap the clumping in tryCatch to skip on error
  tryCatch({
    clumped_chr <- clump_data(
      dat = df_chr,  # you probably want to use df_chr, not bmi3
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 1,
      clump_p2 = 1,
      pop = "EUR",
      bfile = bfile_path,
      plink_bin = genetics.binaRies::get_plink_binary()
    )
    # Append to results only if successful
    clumped_data <- rbind(clumped_data, clumped_chr)
  }, error = function(e) {
    cat("Skipping chromosome", i, "due to error:\n", conditionMessage(e), "\n")
  })
}

exposure<-clumped_data

#Prepare Data for harmonisation

exposure <- exposure %>%
    rename(
      beta.exposure = beta,
      se.exposure = standard_error,
      effect_allele.exposure = effect_allele,
      other_allele.exposure = other_allele,
      eaf.exposure = effect_allele_frequency,
      pval.exposure=p_value,
      N.exposure=n,
      SNP=variant_id
    )
  
  exposure$id.exposure<-"exposure"
  exposure$exposure<-'exposure'
  
  
  
  outcome <- outcome %>%
    rename(
      beta.outcome = beta,
      se.outcome = standard_error,
      effect_allele.outcome = effect_allele,
      other_allele.outcome = other_allele,
      eaf.outcome = effect_allele_frequency,
      pval.outcome=p_value,
      N.outcome=sample_size,
      SNP=variant_id
    )
  
  outcome$id.outcome<-"T1D"
  outcome$outcome<-'T1D'



#HArmonisation
  dat <- harmonise_data(
    exposure_dat = exposure,
    outcome_dat = outcome
  )
#perform proxy search if needed 

# Save Files

saveRDS(exposure,'/n')
saveRDS(dat, '/n')

#Perform MR

results<-mr(dat)
saveRDS(results,'/n')
write.table(results,'/n', col.names=T, row.names=F,quote=F, sep='\t')
