#rsid extraction for Sakaue et al. T1D GWAS

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")

library(usethis)
library(devtools)
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

exposure<-as.data.frame(fread('/n')) #pathway to Sakaue et al. T1D GWAS
exposure <- exposure %>%
    rename(
      CHR = chromosome,
      BP = base_pair_location)
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
data_rsids<-colochelpR::convert_loc_to_rs(df=exposure, dbsnp = dbsnp_144)
