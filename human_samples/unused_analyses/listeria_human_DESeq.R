library(DESeq2)
setwd("/Users/taylorreiter/Desktop/Weimer/Po_Seq/Listeria")

# Read in count data
    counts_all<-read.csv("fc_run_all_humanCOUNTSrenameed.csv")
    filter_out<-c("S10","S9", "S13","S18","S16","S15","S12","S11","PO17", "PO18", "PO19","PO20","PO5","PO6","PO7","PO8","S19", "S14", "S17")
    listeria_counts<-counts_all[, !(names(counts_all) %in% filter_out)]
    head(listeria_counts)
    
# Separate out treatments
    d_treatment6<-c("S7", "S6", "S5", "S8", "S28", "S29", "S30", "S31")
    d_treatment6_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment6)]
    # Round featureCounts output so it will be accepted by DESeq
    d_treatment6_counts<-round(d_treatment6_counts)
    head(d_treatment6_counts)
    # Add colData
    colData6 <- read.csv("colData6.csv", row.names=1)
    # Create dd6 object
    dds6 <- DESeqDataSetFromMatrix(countData = d_treatment6_counts,
                                  colData = colData6,
                                  design = ~ condition)
    # Pre-filter lowly expressed genes
    dds6 <- dds6[ rowSums(counts(dds6)) > 1, ]
    # Preform differential expression analysis (lib size adjustment, dispersion, etc. built in)
    dds6.difex <- DESeq(dds6)
    # Write results
    res6 <- results(dds6.difex)
    # export csv (reorder by padj after export)
    write.csv(res6, file = "deseq_t6.csv")
    
    # Treatment 7
    d_treatment7<-c("S27", "S20", "S21", "S22", "S28", "S29", "S30", "S31")
    d_treatment7_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment7)]
    d_treatment7_counts<-round(d_treatment7_counts)
    head(d_treatment7_counts)
    colData7 <- read.csv("colData7.csv", row.names=1)
    dds7 <- DESeqDataSetFromMatrix(countData = d_treatment7_counts,
                                   colData = colData7,
                                   design = ~ condition)
    dds7 <- dds7[ rowSums(counts(dds7)) > 1, ]
    dds7.difex <- DESeq(dds7)
    res7 <- results(dds7.difex)
    write.csv(res7, file = "deseq_t7.csv")
    
    # Treatment 8
    d_treatment8<-c("S23", "S24", "S25", "S26", "S32", "S33", "S34")
    d_treatment8_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment8)]
    d_treatment8_counts<-round(d_treatment8_counts)
    head(d_treatment8_counts)
    colData8 <- read.csv("colData8.csv", row.names=1)
    dds8 <- DESeqDataSetFromMatrix(countData = d_treatment8_counts,
                                   colData = colData8,
                                   design = ~ condition)
    dds8 <- dds8[ rowSums(counts(dds8)) > 1, ]
    dds8.difex <- DESeq(dds8)
    res8 <- results(dds8.difex)
    write.csv(res8, file = "deseq_t8.csv")
    
    # Treatment 9
    d_treatment9<-c("S27", "S20", "S21", "S22", "S32", "S33", "S34")
    d_treatment9_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment9)]
    d_treatment9_counts<-round(d_treatment9_counts)
    head(d_treatment9_counts)
    colData9 <- read.csv("colData9.csv", row.names=1)
    dds9 <- DESeqDataSetFromMatrix(countData = d_treatment9_counts,
                                   colData = colData9,
                                   design = ~ condition)
    dds9 <- dds9[ rowSums(counts(dds9)) > 1, ]
    dds9.difex <- DESeq(dds9)
    res9 <- results(dds9.difex)
    write.csv(res9, file = "deseq_t9.csv")
    
    # Treatment 10
    d_treatment10<-c("S7", "S6", "S5", "S8", "S23", "S24", "S25", "S26")
    d_treatment10_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment10)]
    d_treatment10_counts<-round(d_treatment10_counts)
    head(d_treatment10_counts)
    colData10 <- read.csv("colData10.csv", row.names=1)
    dds10 <- DESeqDataSetFromMatrix(countData = d_treatment10_counts,
                                   colData = colData10,
                                   design = ~ condition)
    dds10 <- dds10[ rowSums(counts(dds10)) > 1, ]
    dds10.difex <- DESeq(dds10)
    res10 <- results(dds10.difex)
    write.csv(res10, file = "deseq_t10.csv")
    
    # Treatment 11
    d_treatment11<-c("S28", "S29", "S30", "S31", "S32", "S33", "S34")
    d_treatment11_counts<-listeria_counts[, (names(listeria_counts) %in% d_treatment11)]
    d_treatment11_counts<-round(d_treatment11_counts)
    head(d_treatment11_counts)
    colData11 <- read.csv("colData11.csv", row.names=1)
    dds11 <- DESeqDataSetFromMatrix(countData = d_treatment11_counts,
                                   colData = colData11,
                                   design = ~ condition)
    dds11 <- dds11[ rowSums(counts(dds11)) > 1, ]
    dds11.difex <- DESeq(dds11)
    res11 <- results(dds11.difex)
    write.csv(res11, file = "deseq_t11.csv")
