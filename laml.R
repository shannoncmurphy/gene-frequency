# LAML #
LAMLmut.maf <- GDCquery_Maf(tumor = "LAML", pipelines = "mutect2")

laml_clin <- GDCquery_clinic(project = "TCGA-LAML", type = "clinical")

colnames(laml_clin)[1] <- "Tumor_Sample_Barcode"
laml_clin$Overall_Survival_Status <- 1
laml_clin$Overall_Survival_Status[which(laml_clin$vital_status != "dead")] <- 0
laml_clin$time <- laml_clin$days_to_death
laml_clin$time[is.na(laml_clin$days_to_death)] <- laml_clin$days_to_last_follow_up[is.na(laml_clin$days_to_death)]

laml.maf <- read.maf(maf = LAMLmut.maf, clinicalData = laml_clin, isTCGA = T)

