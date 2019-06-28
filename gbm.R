# GBM #
GBMmut.maf <- GDCquery_Maf(tumor = "GBM", pipelines = "mutect2")

gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "clinical")

colnames(gbm_clin)[1] <- "Tumor_Sample_Barcode"
gbm_clin$Overall_Survival_Status <- 1
gbm_clin$Overall_Survival_Status[which(gbm_clin$vital_status != "dead")] <- 0
gbm_clin$time <- gbm_clin$days_to_death
gbm_clin$time[is.na(gbm_clin$days_to_death)] <- gbm_clin$days_to_last_follow_up[is.na(gbm_clin$days_to_death)]

gbm.maf <- read.maf(maf = GBMmut.maf, clinicalData = gbm_clin, isTCGA = T)

gbm_summary <- getGeneSummary(gbm.maf)
gbm_summary 


gbm.gistic <- readGistic(gisticAllLesionsFile = "C://Users//shaan//OneDrive//Documents//Thesis//gbm//all_lesions.conf_99_gbm.txt",
                         gisticAmpGenesFile = "C://Users//shaan//OneDrive//Documents//Thesis//gbm//amp_genes.conf_99_gbm.txt",
                         gisticDelGenesFile = "C://Users//shaan//OneDrive//Documents//Thesis//gbm//del_genes.conf_99_gbm.txt",
                         gisticScoresFile = "C://Users//shaan//OneDrive//Documents//Thesis//gbm//scores_gbm.gistic")

gisticChromPlot(gistic = gbm.gistic,
                markBands = "all")
gisticOncoPlot(gistic = gbm.gistic,
               clinicalData = gbm_clin,
               top = 10)

gbm.genes <- plyr::rbind.fill.matrix(gbm.gistic@gene.summary[["Hugo_Symbol"]],gbm.maf@gene.summary[["Hugo_Symbol"]])

ocurrence <- data.frame(table(gbm.genes))
gbm.genes <- occurance[occurance$Freq > 1,]
gbm.genes[gbm.genes$genes %in% occurance$genes[occurance$Freq >1],]

which(gbm.genes$genes == "ARF1")
which(gbm.maf@gene.summary[["Hugo_Symbol"]] == "RAB11FIP1")
which(gbm.gistic@gene.summary[["Hugo_Symbol"]] == "ARF6")

interest.genes <- c("TNK2","NDRG1","CAV1","IQSEC1","HAX1","EPS15","NUMB","MDM2","DAB2","MYO5A","MYO5B","MYO5C","RAB4A","RAB5A","RAB7A","RAB8A","RAB10","RAB25","RAB35","RAB11FIP1","RAB11FIP2","RAB11FIP3","RAB11FIP4","RAB11FIP15","ARF1","ARF6")

gbm.genes$Freq <- NULL
colnames(gbm.genes) <- "interest.genes"
gbm.genes <- add_row(gbm.genes,interest.genes)
colnames(gbm.genes) <- "Hugo_Symbol"

#mutations
mut.genes <- gbm.maf@gene.summary
gbm.genes$missense_mut <- mut.genes$Missense_Mutation[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$nonsense_mut <- mut.genes$Nonsense_Mutation[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$splice_site <- mut.genes$Splice_Site[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$in_frame_del <- mut.genes$In_Frame_Del[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$in_frame_ins <- mut.genes$In_Frame_Ins[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$frame_shift_ins <- mut.genes$Frame_Shift_Ins[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$frame_shift_del <- mut.genes$Frame_Shift_Del[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
gbm.genes$total_mut <- mut.genes$MutatedSamples[match(gbm.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]

#alterations
alt.genes <- gbm.gistic@gene.summary
gbm.genes$amp <- alt.genes$Amp[match(gbm.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
gbm.genes$del <- alt.genes$Del[match(gbm.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
gbm.genes$total_alt <- alt.genes$AlteredSamples[match(gbm.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]

total_alt_samples <- 577
total_mut_samples <- 392

gbm.genes$alt_freq <- (gbm.genes$total_alt/total_alt_samples)*100
gbm.genes$mut_freq <- (gbm.genes$total_mut/total_mut_samples)*100

gbm.genes$total <- (gbm.genes$alt_freq + gbm.genes$mut_freq)

#max freq
gbm.genes[which.max(gbm.genes$total),]
