#proxy search for Verma et al. T1D GWAS

merge_proxy_with_outcome <- function(ancestry, proxy_MR, outcome, save_path) {
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(LDlinkR)
  
  # Define the population according to the ancestry 
  if (ancestry == "AFR") {
    pop <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
  } else if (ancestry == "ADM") {
    pop <- c("MXL", "PUR", "CLM", "PEL")
  } else if (ancestry == "EUR") {
    pop <- c("CEU", "TSI", "FIN", "GBR","IBS")
  }
  
  
  # Start the proxy search
  token <- 'personalize' 
  LDproxy_batch(proxy_MR$SNP, token = token, genome_build = "grch38", append = TRUE, pop = pop)
  
  # Read the proxy results
  proxySearch <- read.table('combined_query_snp_list_grch38.txt',header=TRUE, sep='\t', row.names=NULL)
  

  outcomeToMerge <- outcome[0, ]
  for (snp_p in proxy_MR$SNP) {
    cat("Processing:", snp_p, "\n")
    
    ldp <- proxySearch[proxySearch$query_snp == snp_p & proxySearch$R2 > 0.8, ]
    ldp <- ldp[ldp$RS_Number %in% outcome$SNP, ]
    
    if (nrow(ldp) == 0) next
    
    ldp2 <- ldp[1, ]
    toMerge <- outcome[outcome$SNP == ldp2$RS_Number, ]
    
    if (ldp2$Distance != 0) {
      alleles <- strsplit(ldp2$Correlated_Alleles, ",")[[1]]
      
      if (length(alleles) < 2) {
        cat("Correlated_Alleles format error for", snp_p, "\n")
        next
      }
      
      ref <- strsplit(alleles[1], "=")[[1]][2]
      alt <- strsplit(alleles[2], "=")[[1]][2]
      
      if (length(toMerge$effect_allele.outcome) == 0 || is.na(ref) || is.na(alt)) {
        cat("Missing allele info for", snp_p, "\n")
        next
      }
      
      # Align alleles
      if (toMerge$effect_allele.outcome == alt) {
        toMerge$effect_allele.outcome <- ref
        toMerge$other_allele.outcome  <- alt
      } else {
        toMerge$effect_allele.outcome <- alt
        toMerge$other_allele.outcome  <- ref
      }
    }
    
    # Match column names and types explicitly
    toMerge <- toMerge[, names(outcomeToMerge), drop = FALSE]
    
    cat("Adding proxy SNP:", toMerge$SNP, "\n")
    
    outcomeToMerge <- rbind(outcomeToMerge, toMerge)
    
  }
  
  # Combine original outcome and proxies
  outcome_final <- rbind(outcome, outcomeToMerge)
  
  # Save
  saveRDS(outcome_final, save_path)
  cat("Outcome with proxies saved to:", save_path, "\n")
  
  return(outcome_final)
}






