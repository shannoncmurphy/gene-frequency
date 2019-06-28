# LAML #
LAMLmut.maf <- GDCquery_Maf(tumor = "LAML", pipelines = "mutect2")

laml_clin <- GDCquery_clinic(project = "TCGA-LAML", type = "clinical")

colnames(laml_clin)[1] <- "Tumor_Sample_Barcode"
laml_clin$Overall_Survival_Status <- 1
laml_clin$Overall_Survival_Status[which(laml_clin$vital_status != "dead")] <- 0
laml_clin$time <- laml_clin$days_to_death
laml_clin$time[is.na(laml_clin$days_to_death)] <- laml_clin$days_to_last_follow_up[is.na(laml_clin$days_to_death)]

laml.maf <- read.maf(maf = LAMLmut.maf, clinicalData = laml_clin, isTCGA = T)

laml_summary <- getGeneSummary(laml.maf)

laml.gistic <- readGistic(gisticAllLesionsFile = "C://Users//shaan//OneDrive//Documents//Thesis//laml//all_lesions.conf_99.txt",
                                        gisticAmpGenesFile = "C://Users//shaan//OneDrive//Documents//Thesis//laml//amp_genes.conf_99.txt",
                                        gisticDelGenesFile = "C://Users//shaan//OneDrive//Documents//Thesis//laml//del_genes.conf_99.txt",
                                        gisticScoresFile = "C://Users//shaan//OneDrive//Documents//Thesis//laml//scores.gistic")

gisticChromPlot(gistic = laml.gistic,
                markBands = "all")

gisticOncoPlot(gistic = laml.gistic,
               clinicalData = laml_clin,
               top = 10)

laml.genes <- plyr::rbind.fill.matrix(laml.gistic@gene.summary[["Hugo_Symbol"]],laml.maf@gene.summary[["Hugo_Symbol"]])

occrance <- data.frame(table(laml.genes))
laml.genes <- occurance[occurance$Freq >1,]
length(laml.genes)
laml.genes$Freq <- NULL

which(laml.genes$genes == "ARF6")

colnames(laml.genes) <- "interest.genes"
laml.genes <- add_row(laml.genes, interest.genes)
colnames(laml.genes) <- "Hugo_Symbol"

mut.genes <- laml.maf@gene.summary
laml.genes$missense_mut <- mut.genes$Missense_Mutation[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$nonsense_mut <- mut.genes$Nonsense_Mutation[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$splice_site <- mut.genes$Splice_Site[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$in_frame_del <- mut.genes$In_Frame_Del[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$in_frame_ins <- mut.genes$In_Frame_Ins[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$frame_shift_del <- mut.genes$Frame_Shift_Del[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
laml.genes$total_mut <- mut.genes$MutatedSamples[match(laml.genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]

alt.genes <- laml.gistic@gene.summary
laml.genes$amp <- alt.genes$Amp[match(laml.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
laml.genes$del <- alt.genes$Del[match(laml.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
laml.genes$total_alt <- alt.genes$AlteredSamples[match(laml.genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]

total_alt_samples <- 191
total_mut_samples <- 140

laml.genes$alt_freq <- (laml.genes$total_alt/total_alt_samples)*100
laml.genes$mut_freq <- (laml.genes$total_mut/total_mut_samples)*100
laml.genes$total <- (laml.genes$alt_freq + laml.genes$mut_freq)

laml.genes[which.max(laml.genes$total),]
