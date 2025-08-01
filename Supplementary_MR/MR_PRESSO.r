

#install MR PRESSO package 

if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")

#import MRPRESSO library 
library(MRPRESSO)

#dat is harmonised data 
#precise the column of your dataframe if changed after haromonisation 
#we're interested in the Outlier-corrected model 
presso_results <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE,
          data = dat)



