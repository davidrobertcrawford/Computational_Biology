# "SA_isoform_TPM_gene_TPM_p_LL_skcm.R"
# 2023-06-11

rm(list=ls())
gc()

library(tidyverse)
library(survival)
library(data.table)
library(foreach)
library(doParallel)

NUM_CORES=12

registerDoParallel(cores=NUM_CORES)
round_size = 24000

source("Scripts/dt_qnorm.R")


# Load gene expression ----------------------------------------------------

load("Data_processed/gencode_v23_dictionary.RData")
pc_genes = unique(gencode[Transcript_type=="protein_coding", Gene_Id])
pc_iso = unique(gencode[Gene_Id %in% pc_genes, Transcript_Id])
gencode_dict = unique(gencode[,.SD,.SDcols=c("Gene_Id", "Transcript_Id")])

gene_dt = fread("Data/Xena/TcgaTargetGtex_rsem_gene_tpm.gz")
gene_pc = gene_dt[sample %in% pc_genes,]
rm(gene_dt)
gc()
gc()

# Load metadata -----------------------------------------------------------

load("Data_processed/TCGA_purity_survival_age.RData")

cancer_vec = "SKCM"

samples_to_keep = unique(purity_survival_age[Cancer_code %in% cancer_vec, Sample_Id])
gene_pc = gene_pc[,.SD,.SDcols=c("sample", samples_to_keep)]
gc()

glimpse(gene_pc[1:3,1:10])
setnames(gene_pc, "sample", "Gene_Id")
gene_melt = melt(gene_pc, id.vars="Gene_Id")
colnames(gene_melt)=c("Gene_Id", "Sample_Id", "TPM")
rm(gene_pc)
gc()

gene_melt[,TPM:=2^(TPM)]
gene_melt[,TPM:=TPM-0.001]
gene_melt[,TPM:=ifelse(TPM<0, 0, TPM)]
setnames(gene_melt, "TPM", "Gene_TPM")


current_cancer="SKCM"

for(current_cancer in cancer_vec){
  
  print(current_cancer)
  cancer_name = unique(purity_survival_age[Cancer_code==current_cancer, Disease])
  load(paste0("Data_processed/Xena/individual_type/TPM_xena_",cancer_name,"_PC.RData"))
  
  samples = purity_survival_age[Cancer_code==current_cancer & sample_type=="Metastatic",]
  sample_ids = unique(samples$Sample_Id)
  
  available_columns = intersect(c("Transcript_Id", sample_ids), colnames(cast_dt))
  
  tx_cancer = cast_dt[Transcript_Id %in% pc_iso,.SD,.SDcols=available_columns]
  
  iso_melt = melt(tx_cancer, id.vars="Transcript_Id")
  colnames(iso_melt)=c("Transcript_Id", "Sample_Id", "TPM")
  setnames(iso_melt, "TPM", "Transcript_TPM")
  iso_melt[,Transcript_Id:=as.character(Transcript_Id)]
  iso_melt = merge(iso_melt, gencode_dict, by="Transcript_Id")
  
  gene_short = gene_melt[Sample_Id %in% sample_ids,]
  
  iso_gene_melt = merge(iso_melt, gene_short, by=c("Sample_Id", "Gene_Id"))
  
  
  meta_num = purity_survival_age[Cancer_code==current_cancer, .SD,.SDcols=c("Sample_Id", "ESTIMATE", "Age")]
  meta_cat = purity_survival_age[Cancer_code==current_cancer, .SD,.SDcols=c("Sample_Id", "EVENT", "TIME_TO_EVENT", "SEX", "sample_type")]
  gene_vec = unique(iso_gene_melt$Transcript_Id)
  
  #gene_vec = gene_vec[1:1000]
  
  xena_gene = merge(iso_gene_melt, meta_num, by="Sample_Id")
  setkey(xena_gene, Transcript_Id)
  rm(gene_short, iso_melt, cast_dt, iso_gene_melt)
  gc()
  gc()
  
  # Get results for each gene --------------------------------------------
  matrix_dt = data.table(Transcript_Id = character(),
                         Cancer_Type=character(),
                         coxph_coef_Transcript_TPM=numeric(),
                         coxph_Pr_Transcript_TPM=numeric(),
                         Gene_Id = character(),
                         coxph_coef_Gene_TPM=numeric(),
                         coxph_Pr_Gene_TPM=numeric(),
                         Num_events = integer(),
                         coxph_coef_Purity=numeric(),
                         coxph_Pr_Purity=numeric(),
                         coxph_coef_Age=numeric(),
                         coxph_Pr_Age=numeric(),
                         logLik_initial=numeric(),
                         logLik_final=numeric(),
                         df=numeric())
  
  num_gene = length(gene_vec)
  gene_rounds = as.integer(num_gene/round_size)
  gene_last = num_gene%%round_size

  for(round in 1:(gene_rounds+1)){
    print(paste0("......starting round ", round))
    start = 1 + (round-1)*round_size
    if(round <= gene_rounds) end=round*round_size
    if(round > gene_rounds) end = start + gene_last - 1
    
    gene_sub_chunks = foreach(gene = gene_vec[start:end], .combine = rbind, 
                              .multicombine = TRUE, .packages="data.table") %dopar%{
                                gene_dt = xena_gene[Transcript_Id==gene,]
                                gene_dt[,Transcript_Id:=NULL]
                                gene_dt[,Gene_Id:=NULL]
                                gnorm_dt = dt_qnorm(gene_dt, "Sample_Id")
                                gnorm_dt[,Transcript_Id:=gene]
                                ranked_SA = merge(gnorm_dt, meta_cat, 
                                                  by="Sample_Id")
                                temp_output = try(coxph(Surv(TIME_TO_EVENT, 
                                                             EVENT) ~ Transcript_TPM +
                                                          Gene_TPM +
                                                          ESTIMATE +
                                                          Age+
                                                          strata(SEX),
                                                        data=ranked_SA))
                                if (class(temp_output) == "try-error"){
                                  temp_row = vector("list", 15)
                                } else if (temp_output$nevent==0){
                                  temp_row = vector("list", 15)
                                }else {
                                  tsum = summary(temp_output)
                        
                                  temp_output2 = try(coxph(Surv(TIME_TO_EVENT, 
                                                                EVENT) ~ 
                                                             #Gene_TPM +
                                                             ESTIMATE +
                                                             Age+
                                                             strata(SEX),
                                                           data=ranked_SA))
                                  if (class(temp_output2) == "try-error"){
                                    temp_row = vector("list", 15)
                                  } else if (temp_output2$nevent==0){
                                    temp_row = vector("list", 15)
                                  }else {
                                    tsum2 = summary(temp_output2)
                                    temp_row = list(gene, current_cancer,
                                                    tsum[['coefficients']][1],
                                                    tsum[['coefficients']][17],
                                                    gencode_dict[Transcript_Id==gene, Gene_Id],
                                                    tsum[['coefficients']][2],
                                                    tsum[['coefficients']][18],
                                                    tsum$nevent,
                                                    tsum[['coefficients']][3],
                                                    tsum[['coefficients']][19],
                                                    tsum[['coefficients']][4],
                                                    tsum[['coefficients']][20],
                                                    tsum2[['loglik']][[2]],
                                                    tsum[['loglik']][[2]],
                                                    2
                                    )
                                  }
                                }
                                temp_row
                              }
    mat1 = as.data.table(gene_sub_chunks)
    premat = lapply(mat1, unlist)
    matrix_temp = as.data.table(premat)
    colnames(matrix_temp) = c("Transcript_Id", "Cancer_Type", "coxph_coef_Transcript_TPM",
                              "coxph_Pr_Transcript_TPM", "Gene_Id", "coxph_coef_Gene_TPM",
                              "coxph_Pr_Gene_TPM", "Num_events",
                              "coxph_coef_Purity", "coxph_Pr_Purity",
                              "coxph_coef_Age", "coxph_Pr_Age", 
                              "logLik_initial", "logLik_final", "df")
    matrix_temp = matrix_temp[!is.na(matrix_temp$coxph_coef_Transcript_TPM),]
    matrix_dt = rbind(matrix_dt, matrix_temp)
    rm(matrix_temp, mat1, premat)
    gc()
  }
  
  
  dir_phrase="SA_iso_gene_purity_LL"
  
  save(matrix_dt, file=paste0("Results/SA_skcm/",dir_phrase,"_",current_cancer,".RData"))

}
