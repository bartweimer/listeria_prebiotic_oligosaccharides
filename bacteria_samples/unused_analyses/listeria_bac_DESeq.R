setwd("/Users/taylorreiter/Desktop/Weimer/Po_Seq/Listeria/bac_listeria")

library(DESeq2)
# Read in count data
counts_bac<-read.csv("fc_run_all_listeriaCOUNTSrenamed.csv", row.name=1)
head(counts_bac)


# Separate out treatments
      list_treatment6<-c("S7", "S6", "S5", "S8", "S28", "S29", "S30", "S31")
      list_treatment6_counts<-counts_bac[, (names(counts_bac) %in% list_treatment6)]
      # Round featureCounts output so it will be accepted by DESeq
      list_treatment6_counts<-round(list_treatment6_counts)
      head(list_treatment6_counts)
      # Add colData
      list_colData6 <- read.csv("list_colData6.csv", row.names=1)
      # Create dd6 object
      list_dds6 <- DESeqDataSetFromMatrix(countData = list_treatment6_counts,
                                     colData = list_colData6,
                                     design = ~ condition)
      # Preform differential expression analysis (lib size adjustment, dispersion, etc. built in)
      list_dds6.difex <- DESeq(list_dds6)
      # Write results
      list_res6 <- results(list_dds6.difex)
      # export csv (reorder by padj after export)
      write.csv(list_res6, file = "deseq_list_t6.csv")

      # Treatment 7
      list_treatment7<-c("S27", "S20", "S21", "S22", "S28", "S29", "S30", "S31")
      list_treatment7_counts<-counts_bac[, (names(counts_bac) %in% list_treatment7)]
      list_treatment7_counts<-round(list_treatment7_counts)
      head(list_treatment7_counts)
      list_colData7 <- read.csv("list_colData7.csv", row.names=1)
      list_dds7 <- DESeqDataSetFromMatrix(countData = list_treatment7_counts,
                                     colData = list_colData7,
                                     design = ~ condition)
      list_dds7.difex <- DESeq(list_dds7)
      list_res7 <- results(list_dds7.difex)
      write.csv(list_res7, file = "deseq_list_t7.csv")
      
      # Treatment 8
      list_treatment8<-c("S23", "S24", "S25", "S26", "S32", "S33", "S34")
      list_treatment8_counts<-counts_bac[, (names(counts_bac) %in% list_treatment8)]
      list_treatment8_counts<-round(list_treatment8_counts)
      head(list_treatment8_counts)
      list_colData8 <- read.csv("list_colData8.csv", row.names=1)
      list_dds8 <- DESeqDataSetFromMatrix(countData = list_treatment8_counts,
                                     colData = list_colData8,
                                     design = ~ condition)
      list_dds8.difex <- DESeq(list_dds8)
      list_res8 <- results(list_dds8.difex)
      write.csv(list_res8, file = "deseq_list_t8.csv")
      
      # Treatment 9
      list_treatment9<-c("S27", "S20", "S21", "S22", "S32", "S33", "S34")
      list_treatment9_counts<-counts_bac[, (names(counts_bac) %in% list_treatment9)]
      list_treatment9_counts<-round(list_treatment9_counts)
      head(list_treatment9_counts)
      list_colData9 <- read.csv("list_colData9.csv", row.names=1)
      list_dds9 <- DESeqDataSetFromMatrix(countData = list_treatment9_counts,
                                     colData = list_colData9,
                                     design = ~ condition)
      list_dds9.difex <- DESeq(list_dds9)
      list_res9 <- results(list_dds9.difex)
      write.csv(list_res9, file = "deseq_list_t9.csv")
      
      # Treatment 10
      list_treatment10<-c("S7", "S6", "S5", "S8", "S23", "S24", "S25", "S26")
      list_treatment10_counts<-counts_bac[, (names(counts_bac) %in% list_treatment10)]
      list_treatment10_counts<-round(list_treatment10_counts)
      head(list_treatment10_counts)
      list_colData10 <- read.csv("list_colData10.csv", row.names=1)
      list_dds10 <- DESeqDataSetFromMatrix(countData = list_treatment10_counts,
                                      colData = list_colData10,
                                      design = ~ condition)
      list_dds10.difex <- DESeq(list_dds10)
      list_res10 <- results(list_dds10.difex)
      write.csv(list_res10, file = "deseq_list_t10.csv")

      # Treatment 11
      list_treatment11<-c("S28", "S29", "S30", "S31", "S32", "S33", "S34")
      list_treatment11_counts<-counts_bac[, (names(counts_bac) %in% list_treatment11)]
      list_treatment11_counts<-round(list_treatment11_counts)
      head(list_treatment11_counts)
      list_colData11 <- read.csv("list_colData11.csv", row.names=1)
      list_dds11 <- DESeqDataSetFromMatrix(countData = list_treatment11_counts,
                                      colData = list_colData11,
                                      design = ~ condition)
      list_dds11.difex <- DESeq(list_dds11)
      list_res11 <- results(list_dds11.difex)
      write.csv(list_res11, file = "deseq_list_t11.csv")
      