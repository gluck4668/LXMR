
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXMR")
library(LXMR)

rm(list=ls())





#----local data TwosampleMR---

#-----inpute exposure data---
{exposure="Depression"
  exp_dat="finngen_R10_F5_DEPRESSION_RECURRENT.gz"
  exp_p=5e-8

  snp_exp = "rsids"
  beta_exp = "beta"
  se_exp = "sebeta"
  effect_allele_exp = "alt"
  other_allele_exp = "ref"
  eaf_exp = "af_alt"
  pval_exp = "pval"

  clump_kb = 10000
  clump_r2 = 0.001
}

#----inpute outcome data ---
{outcome="NAFLD"
  out_dat="finngen_R10_NAFLD.gz"

  snp_out = "rsids"
  beta_out = "beta"
  se_out = "sebeta"
  effect_allele_out = "alt"
  other_allele_out = "ref"
  eaf_out = "af_alt"
  pval_out = "pval"}

#---run the main function---

devtools::load_all()

LXMR(exposure,exp_dat,exp_p,snp_exp,beta_exp,se_exp,
     effect_allele_exp,other_allele_exp,eaf_exp,pval_exp,
     clump_kb,clump_r2,outcome,out_dat,snp_out,beta_out,
     se_out,effect_allele_out,other_allele_out,eaf_out,pval_out)
