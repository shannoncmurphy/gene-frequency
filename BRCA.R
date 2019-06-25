library(maftools)
library(TCGAbiolinks)
library(tibble)

BRCAmut.maf <- GDCquery_Maf(tumor = "BRCA", pipelines = "mutect2")

brca_clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "Clinical")

colnames(brca_clin)[1] <- "Tumor_Sample_Barcode"
brca_clin$Overall_Survival_Status <- 1
brca_clin$Overall_Survival_Status[which(brca_clin$vital_status != "dead")] <- 0
brca_clin$time <- brca_clin$days_to_death
brca_clin$time[is.na(brca_clin$days_to_death)] <- brca_clin$days_to_last_follow_up[is.na(brca_clin$days_to_death)]

brca.maf <- read.maf(maf = BRCAmut.maf, clinicalData = brca_clin, isTCGA = T)

brca_summary <- getGeneSummary(brca.maf)

all.lesions <- system.file("extdata","all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata","amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata","scores.gistic", package = "maftools")
brca.gistic <- readGistic(gisticAllLesionsFile = all.lesions,
                          gisticAmpGenesFile = amp.genes,
                          gisticDelGenesFile = del.genes,
                          gisticScoresFile = scores.gis)

brca.gistic

gisticChromPlot(gistic = brca.gistic,
                markBands = "all")

gisticOncoPlot(gistic = brca.gistic,
               clinicalData = brca_clin,
               top = 10)

#list of all genes
genes <- plyr::rbind.fill.matrix(brca.gistic@gene.summary[["Hugo_Symbol"]],brca.maf@gene.summary[["Hugo_Symbol"]])

#only keep genes that occur in both
occurance <- data.frame(table(genes))
genes <- occurance[occurance$Freq > 1,]
genes[genes$genes %in% occurance$genes[occurance$Freq > 1],]

#remove freq column
genes$Freq <- NULL

which(genes$genes == "ARF6")
which(brca.maf@gene.summary[["Hugo_Symbol"]] == "ARF6")
which(brca.gistic@gene.summary[["Hugo_Symbol"]] == "RAB35")

#genes of interest
interest.genes <- c("TNK2","CAV1","IQSEC1","HAX1","EPS15","NUMB","MDM2","DAB2","MYO5A","MYO5B","MYO5C","RAB4A","RAB5A","RAB7A","RAB8A","RAB10","RAB25","RAB35","RAB11FIP1","RAB11FIP2","RAB11FIP3","RAB11FIP4","RAB11FIP15","ARF1","ARF6")

#Add genes of interest
colnames(genes) <- "interest.genes"
genes <- add_row(genes,interest.genes)
colnames(genes) <- "Hugo_Symbol"

#Add mutation data matching to genes
mut.genes <- brca.maf@gene.summary
genes$missense_mut <- mut.genes$Missense_Mutation[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$nonsense_mut <- mut.genes$Nonsense_Mutation[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$splice_site <- mut.genes$Splice_Site[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$in_frame_del <- mut.genes$In_Frame_Del[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$in_frame_ins <- mut.genes$In_Frame_Ins[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$frame_shift_ins <- mut.genes$Frame_Shift_Ins[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$frame_shift_del <- mut.genes$Frame_Shift_Del[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]
genes$total_mut <- mut.genes$MutatedSamples[match(genes$Hugo_Symbol,mut.genes$Hugo_Symbol)]

#Add alteration data matching to genes
alt.genes <- brca.gistic@gene.summary
genes$amp <- alt.genes$Amp[match(genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
genes$del <- alt.genes$Del[match(genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]
genes$total_alt <- alt.genes$AlteredSamples[match(genes$Hugo_Symbol,alt.genes$Hugo_Symbol)]

total_alt_samples <- 1080
total_mut_samples <- 985

#find frequency of mutations and alterations
genes$alt_freq <- (genes$total_alt/total_alt_samples)*100
genes$mut_freq <- (genes$total_mut/total_mut_samples)*100

#total frequency
genes$total <- (genes$alt_freq + genes$mut_freq)

#find gene with highest frequency
genes[which.max(genes$total),]
