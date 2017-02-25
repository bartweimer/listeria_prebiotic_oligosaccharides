# Set working directory (location of relevant files)
    setwd("/Users/taylorreiter/Desktop/Weimer/Po_seq/Listeria_human_counts_output")

# load required packages
    library(edgeR)
    library(limma)

    setwd("/Users/taylorreiter/Desktop/Weimer/Po_Seq/Listeria/human_listeria/")
# Import count tables. 
# Import total counts for all experiments and subset to relevant later
    counts_all<-read.table("fc_run_all_humanCOUNTS.txt")

# Change column names from file name to sample name
    previous.names<-colnames(counts_all)
    samples_names<-c("S10","S7","S6","S5","S9","S23","S13","S24","S25","S18","S16","S26","S27","S20","S21","S22","S28","S29","S15","S30","S12","S31","S11","S32", "PO17", "PO18", "PO19","PO20","PO5","PO6","PO7","PO8","S33", "S17", "S8", "S19", "S14", "S34")
    counts_all<-as.data.frame(counts_all)
    library(data.table)
    setnames(counts_all, previous.names, samples_names)

# Export as csv so renaming won't need to be repeated
    write.csv(counts_all, file="fc_run_all_humanCOUNTSrenameed.csv")

# Subset to listeria
    filter_out<-c("S10","S9", "S13","S18","S16","S15","S12","S11","PO17", "PO18", "PO19","PO20","PO5","PO6","PO7","PO8","S19", "S14", "S17")
    listeria_counts<-counts_all[, !(names(counts_all) %in% filter_out)]

# Subset for each treatment & create a dataframe; create group object for treatment
    # Treatment 6
    treatment6<-c("S7", "S6", "S5", "S8", "S28", "S29", "S30", "S31")
    treatment6_counts<-listeria_counts[, (names(listeria_counts) %in% treatment6)]
    head(treatment6_counts)
    group6<-c("control", "control", "control", "listeria", "listeria", "listeria", "listeria", "control")
    # Treatment 7
    treatment7<-c("S27", "S20", "S21", "S22", "S28", "S29", "S30", "S31")
    treatment7_counts<-listeria_counts[, (names(listeria_counts) %in% treatment7)]
    head(treatment7_counts)
    group7<-c("control", "control", "control", "control", "listeria", "listeria", "listeria", "listeria")
    # Treatment 7 sans S20
    treatment7<-c("S27", "S20", "S21", "S22", "S28", "S29", "S30", "S31")
    treatment7_counts<-listeria_counts[, (names(listeria_counts) %in% treatment7)]
    head(treatment7_counts)
    group7<-c("control", "control", "control", "control", "listeria", "listeria", "listeria", "listeria")
    # Treatment 8
    treatment8<-c("S23", "S24", "S25", "S26", "S32", "S33", "S34")
    treatment8_counts<-listeria_counts[, (names(listeria_counts) %in% treatment8)]
    head(treatment8_counts)
    group8<-c("control", "control", "control", "control", "listeria", "listeria", "listeria")
    # Treatment 9
    treatment9<-c("S27", "S20", "S21", "S22", "S32", "S33", "S34")
    treatment9_counts<-listeria_counts[, (names(listeria_counts) %in% treatment9)]
    head(treatment9_counts)
    group9<-c("control", "control", "control", "control", "listeria", "listeria", "listeria")
    # Treatment 10
    treatment10<-c("S7", "S6", "S5", "S8", "S23", "S24", "S25", "S26")
    treatment10_counts<-listeria_counts[, (names(listeria_counts) %in% treatment10)]
    head(treatment10_counts)
    group10<-c("control", "control", "control", "listeria", "listeria", "listeria", "listeria", "control")
    # Treatment 11
    treatment11<-c("S28", "S29", "S30", "S31", "S32", "S33", "S34")
    treatment11_counts<-listeria_counts[, (names(listeria_counts) %in% treatment11)]
    head(treatment11_counts)
    group11<-c("control", "control", "control", "control", "listeria", "listeria", "listeria")

# Add library size vectors for normalization
    # make objects for size of each library
    S7.lib<-14802076
    S6.lib<-25476648
    S5.lib<-35610664 
    S8.lib<-52401510 
    S23.lib<-102661582 
    S24.lib<-18693011
    S25.lib<-63459929 
    S26.lib<-16689093  
    S27.lib<-18567572  
    S20.lib<-69055612  
    S21.lib<-67387756  
    S22.lib<-34129549  
    S28.lib<-67209403  
    S29.lib<-79197902  
    S30.lib<-57654165 
    S31.lib<-24538213 
    S32.lib<-51888906  
    S33.lib<-62610049 
    S34.lib<-28880047 
    
# Make vectors of library size for each treament
    library6<-c(S7.lib, S6.lib, S5.lib, S28.lib, S29.lib, S30.lib, S31.lib, S8.lib)
    library7<-c(S27.lib, S20.lib, S21.lib, S22.lib, S28.lib, S29.lib, S30.lib, S31.lib)
    library8<-c(S23.lib, S24.lib, S25.lib, S26.lib, S32.lib, S33.lib, S34.lib)
    library9<-c(S27.lib, S20.lib, S21.lib, S22.lib, S32.lib, S33.lib, S34.lib)
    library10<-c(S7.lib, S6.lib, S5.lib, S23.lib, S24.lib, S25.lib, S26.lib, S8.lib)
    library11<-c(S28.lib, S29.lib, S30.lib, S31.lib, S32.lib, S33.lib, S34.lib)
    
# Make DGEList with library size
    treatment6_counts<-DGEList(counts=treatment6_counts, group = group6, lib.size = library6)
    treatment7_counts<-DGEList(counts=treatment7_counts, group = group7, lib.size = library7)
    treatment8_counts<-DGEList(counts=treatment8_counts, group = group8, lib.size = library8)
    treatment9_counts<-DGEList(counts=treatment9_counts, group = group9, lib.size = library9)
    treatment10_counts<-DGEList(counts=treatment10_counts, group = group10, lib.size = library10)
    treatment11_counts<-DGEList(counts=treatment11_counts, group = group11, lib.size = library11)

# Filter out lowly expressed genes
# Criterion: normalized to counts per million, 6 or 7 reads per sample in 2 or more libraries per group
    keep6<-rowSums(cpm(treatment6_counts)>1) >=2
    treatment6_counts <- treatment6_counts[keep6, , keep.lib.sizes=FALSE]
    keep7<-rowSums(cpm(treatment7_counts)>1) >=2
    treatment7_counts <- treatment7_counts[keep7, , keep.lib.sizes=FALSE]
    keep8<-rowSums(cpm(treatment8_counts)>1) >=2
    treatment8_counts <- treatment8_counts[keep8, , keep.lib.sizes=FALSE]
    keep9<-rowSums(cpm(treatment9_counts)>1) >=2
    treatment9_counts <- treatment9_counts[keep9, , keep.lib.sizes=FALSE]
    keep10<-rowSums(cpm(treatment10_counts)>1) >=2
    treatment10_counts <- treatment10_counts[keep10, , keep.lib.sizes=FALSE]
    keep11<-rowSums(cpm(treatment11_counts)>1) >=2
    treatment11_counts <- treatment11_counts[keep11, , keep.lib.sizes=FALSE]

# Normalize reads
# Normalized to library size
    treatment6_counts <- calcNormFactors(treatment6_counts)
    treatment7_counts <- calcNormFactors(treatment7_counts)
    treatment8_counts <- calcNormFactors(treatment8_counts)
    treatment9_counts <- calcNormFactors(treatment9_counts)
    treatment10_counts <- calcNormFactors(treatment10_counts)
    treatment11_counts <- calcNormFactors(treatment11_counts)
# Create design vector
    design6<-model.matrix(~group6)
    design7<-model.matrix(~group7)
    design8<-model.matrix(~group8)
    design9<-model.matrix(~group9)
    design10<-model.matrix(~group10)
    design11<-model.matrix(~group11)

# Estimate dispersion
    library(locfit) # dependency not installed when edgeR was installed
    treatment6_counts<-estimateDisp(treatment6_counts, design6)
    treatment7_counts<-estimateDisp(treatment7_counts, design7)
    treatment8_counts<-estimateDisp(treatment8_counts, design8)
    treatment9_counts<-estimateDisp(treatment9_counts, design9)
    treatment10_counts<-estimateDisp(treatment10_counts, design10)
    treatment11_counts<-estimateDisp(treatment11_counts, design11)

# perform exact test for differential expression
# use exactTest() for pairwise comparisons. Based on qCML methods. Strong parallels to Fisher's exact test. 
# this test is only applicable to experiments with a single factor. 
    et6 <- exactTest(treatment6_counts, pair=c("control", "listeria"))
    tags6<-topTags(et6, sort.by="PValue", p.value=1, n=30000)
    write.csv(tags6[[1]], file="edgeR_human_list_p1_t6.csv")
    
    et7 <- exactTest(treatment7_counts, pair=c("control", "listeria"))
    tags7<-topTags(et7, sort.by="PValue", p.value=1, n=30000)
    write.csv(tags7[[1]], file="edgeR_human_list_p1_t7.csv")
    
    et8 <- exactTest(treatment8_counts, pair=c("control", "listeria"))
    tags8<-topTags(et8, sort.by="PValue", p.value=.05, n=30000)
    write.csv(tags8[[1]], file="edgeR_human_list_p1_t8.csv")
    
    et9 <- exactTest(treatment9_counts, pair=c("control", "listeria"))
    tags9<-topTags(et9, sort.by="PValue", p.value=.05, n=30000)
    write.csv(tags9[[1]], file="edgeR_human_list_p1_t9.csv")
    
    et10 <- exactTest(treatment10_counts, pair=c("control", "listeria"))
    tags10<-topTags(et10, sort.by="PValue", p.value=.05, n=30000)
    write.csv(tags10[[1]], file="edgeR_human_list_p1_t10.csv")
    
    et11 <- exactTest(treatment11_counts, pair=c("control", "listeria"))
    tags11<-topTags(et11, sort.by="PValue", p.value=.05, n=30000)
    write.csv(tags11[[1]], file="edgeR_human_list_p1_t11.csv")
   
    plotMDS(listeria_counts)
--------------------------------------------------------
# limma differential expression
# use edgeR inputs, and feed into limma
library(limma)
    treatment6_logCPM <- cpm(treatment6_counts, log=TRUE, prior.count=3)
    fit6 <- lmFit(treatment6_logCPM, design6)
    fit6 <- eBayes(fit6, trend=TRUE)
    table6<-topTable(fit6, coef=ncol(design6), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table6, file="limma_difex_t6.csv")
      
    treatment7_logCPM <- cpm(treatment7_counts, log=TRUE, prior.count=3)
    fit7 <- lmFit(treatment7_logCPM, design7)
    fit7 <- eBayes(fit7, trend=TRUE)
    table7<-topTable(fit7, coef=ncol(design7), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table7, file="limma_difex_t7.csv")
    
    treatment8_logCPM <- cpm(treatment8_counts, log=TRUE, prior.count=3)
    fit8 <- lmFit(treatment8_logCPM, design8)
    fit8 <- eBayes(fit8, trend=TRUE)
    table8<-topTable(fit8, coef=ncol(design8), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table8, file="limma_difex_t8.csv")
    
    treatment9_logCPM <- cpm(treatment9_counts, log=TRUE, prior.count=3)
    fit9 <- lmFit(treatment9_logCPM, design9)
    fit9 <- eBayes(fit9, trend=TRUE)
    table9<-topTable(fit9, coef=ncol(design9), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table9, file="limma_difex_t9.csv")
    
    treatment10_logCPM <- cpm(treatment10_counts, log=TRUE, prior.count=3)
    fit10 <- lmFit(treatment10_logCPM, design10)
    fit10 <- eBayes(fit10, trend=TRUE)
    table10<-topTable(fit10, coef=ncol(design10), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table10, file="limma_difex_t10.csv")
    
    treatment11_logCPM <- cpm(treatment11_counts, log=TRUE, prior.count=3) 
    fit11 <- lmFit(treatment11_logCPM, design11)
    fit11 <- eBayes(fit11, trend=TRUE)
    table11<-topTable(fit11, coef=ncol(design11), sort.by="P", adjust.method="BH", p.value=.05, number=10000)
    write.csv(table11, file="limma_difex_t11.csv")
    
  