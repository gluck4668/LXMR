
local_data <- function(exposure,exp_dat,exp_p,snp_exp,beta_exp,se_exp,
                       effect_allele_exp,other_allele_exp,eaf_exp,pval_exp,
                       clump_kb,clump_r2,outcome,out_dat,snp_out,beta_out,
                       se_out,effect_allele_out,other_allele_out,eaf_out,pval_out){

  #-------R packages ----
  inst_pack <- function(){

    packs <- c("devtools","remote","openxlsx","dplyr","readr","purrr","ggplot2",
               "tidyr","forestploter","tidyverse","stringr","shiny","remotes",
               "MendelianRandomization","simex","data.table","VariantAnnotation")

    pack_uninst <- packs[!packs %in% installed.packages()[,1]]

    if(length(pack_uninst)>0){
      tryCatch(install.packages(pack_unist),error=function(e){e})
      tryCatch(BiocManager::install(pack_unist),error=function(e){e})
    }

    mr_packs <- c("MRInstruments","TwoSampleMR","gwasvcf","gwasglue",
                  "RMVMR","MRPRESSO")
    # "MRInstruments","TwoSampleMR","gwasvcf","gwasglue" belong to MRCIEU;
    # "RMVMR" belong to WSpiller;
    # "MRPRESSO" belong to "rondolab/MR-PRESSO"

    mr_uninst <- mr_packs[!mr_packs %in% installed.packages()[,1]]

    if(length(mr_uninst)>0){
      for(i in mr_packs){
        print(i)
        tryCatch(remotes::install_github(paste0("MRCIEU/",i)),error=function(e){e})
        tryCatch(remotes::install_github(paste0("WSpiller/",i)),error=function(e){e})}
      tryCatch(remotes::install_github("rondolab/MR-PRESSO"),error=function(e){e})
    }

    for(i in c(packs,mr_packs)){
      library(i,character.only = T)}

  }

  inst_pack()

  #---create dir---
  if(!dir.exists("analysis results"))
    dir.create("analysis results")

  #---- Exposure data ---
  exp_file_type <-str_extract(exp_dat,"(?<=\\.).*")
  if(tolower(exp_file_type)=="gz")
    exp_df <- fread(exp_dat) else
    {if(tolower(exp_file_type)=="txt")
      exp_df <- read_table(exp_dat) else
        exp_df <- eval(str2expression(paste0("read.",exp_file_type,"(exp_dat)")))
    }

  exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")")))
  exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)")))
  if(nrow(exp_df)==0){stop(paste("Unpon the condition of pval<",exp_p,", there is no exposure data."))}
  exp_file <- paste0("analysis results/",stringr::str_extract(exp_dat,".*(?=[.])"),".csv")
  write.csv(exp_df,exp_file)

  exp_data <- read_exposure_data(filename = exp_file,
                                 sep = ",",
                                 snp_col = snp_exp,
                                 beta_col = beta_exp,
                                 se_col = se_exp,
                                 effect_allele_col = effect_allele_exp,
                                 other_allele_col = other_allele_exp,
                                 eaf_col = eaf_exp,
                                 pval_col = pval_exp
  )

  exp_data <- TwoSampleMR::clump_data(exp_data,
                                      clump_kb = clump_kb,
                                      clump_r2 = clump_r2)

  exp_snp <- data.frame(exp_data$SNP)
  names(exp_snp) <- snp_out

  #---- Outcome data ---
  out_file_type <-str_extract(out_dat,"(?<=\\.).*")
  if(tolower(out_file_type)=="gz")
    out_df <- fread(out_dat) else
    {if(tolower(out_file_type)=="txt")
      out_df <- read_table(out_dat) else
        out_df <- eval(str2expression(paste0("read.",out_file_type,"(out_dat)")))
    }

  out_df <- merge(exp_snp,out_df,by=snp_out)
  out_df <- eval(str2expression(paste0("distinct(.data = out_df,",snp_out,",.keep_all = T)")))
  out_file <- paste0("analysis results/",stringr::str_extract(out_dat,".*(?=[.])"),".csv")
  write.csv(out_df,out_file)

  out_data <- read_outcome_data(
    snps=exp_data$SNP,
    filename=out_file,
    sep = ",",
    snp_col = snp_out,
    beta_col = beta_out,
    se_col = se_out,
    effect_allele_col = effect_allele_out,
    other_allele_col = other_allele_out,
    eaf_col = eaf_out,
    pval_col = pval_out
  )

  local_data <- list(exp_data=exp_data,out_data=out_data)

  return(local_data)

}
