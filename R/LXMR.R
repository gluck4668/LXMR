
LXMR <- function(exposure,exp_dat,exp_p,snp_exp,beta_exp,se_exp,
                 effect_allele_exp,other_allele_exp,eaf_exp,pval_exp,
                 clump_kb,clump_r2,outcome,out_dat,snp_out,beta_out,
                 se_out,effect_allele_out,other_allele_out,eaf_out,pval_out){

local_data <- local_data(exposure,exp_dat,exp_p,snp_exp,beta_exp,se_exp,
                         effect_allele_exp,other_allele_exp,eaf_exp,pval_exp,
                         clump_kb,clump_r2,outcome,out_dat,snp_out,beta_out,
                         se_out,effect_allele_out,other_allele_out,eaf_out,pval_out)

exp_data <- local_data$exp_data
out_data <- local_data$out_data

mr_analysis <- mr_analysis(exp_data,out_data)

print("The results can be found in the folder of <analysis results> ")

mr_analysis$pict

}
