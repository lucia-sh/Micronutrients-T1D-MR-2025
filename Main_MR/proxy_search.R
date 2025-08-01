Library(LDlinkR)
token="####"


#get the missing SNPs from outcome

proxy_MR<-exposure[!exposure$SNP %in% outcome$SNP,]

# do the proxy search using LDLinkR, a file gonna be generated

LDproxy_batch(proxy_MR$SNP, token = token, genome_build = "grch38", append = TRUE)


#read the generated file 
proxySearch <- read.delim('combined_query_snp_list_grch38.txt', row.names = NULL)

#map the missing SNPs

outcomeToMerge <- NULL
for(snp_p in proxy_MR$SNP){
  print(snp_p)
  
  
  ldp=proxySearch[proxySearch$query_snp==snp_p,]
  ldp <- ldp[ldp$R2 > 0.8,]
  ldp=ldp[ldp$RS_Number%in%outcome$variant_id,]
  if(nrow(ldp)==0){
    next
  }
  ldp2=ldp[1,]
  toMerge=outcome[outcome$variant_id==ldp2$RS_Number,]
  if(ldp2$Distance!=0){
    ref=colsplit(ldp2$Correlated_Alleles,',',c('a','b'))
    alt=colsplit(ref$b,'=',c('exp','out'))
    alt[alt==TRUE]='T'
    ref=colsplit(ref$a,'=',c('exp','out'))
    ref[ref==TRUE]='T'
    # Check that both allele values exist before comparing
    if(length(toMerge$EA) == 0 || length(alt$out) == 0){
      cat("Allele information missing for SNP", snp_p, "\n")
      next
    }
    if(toMerge$EA==alt$out){
      toMerge$EA=alt$exp
      toMerge$Allele1=ref$exp
    }else{
      toMerge$EA=ref$exp
      toMerge$Allele1=alt$exp
    }
  }
  print("Binding...")
  outcomeToMerge=rbind(outcomeToMerge,toMerge)
  saveRDS(outcomeToMerge,"outcomeToMerge.rds")
}

outcome_updated <- rbind(outcome, outcomeToMerge)
outcome <- outcome_updated