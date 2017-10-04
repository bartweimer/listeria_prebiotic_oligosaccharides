setwd("/Users/bcweimer/Desktop/Po_seq/hisat_output/")
library(Rsubread)

fls <- dir(".", "bam$")
fls


fc_human<-featureCounts(files=fls,
                           # annotation
                           annot.inbuilt=NULL,
                           annot.ext="Homo_sapiens.GRCh38.86.gtf",
                           isGTFAnnotationFile=TRUE,
                           GTF.featureType="gene",
                           GTF.attrType="gene_name",
                           chrAliases=NULL,
                           
                           # level of summarization
                           useMetaFeatures=FALSE,
                           
                           # overlap between reads and features
                           allowMultiOverlap=TRUE,
                           minOverlap=10,
                           largestOverlap=FALSE,
                           readExtension5=0,
                           readExtension3=0,
                           read2pos=NULL,
                           
                           # multi-mapping reads
                           countMultiMappingReads=TRUE,
                           fraction=TRUE,
                           
                           # read filtering
                           minMQS=0,
                           splitOnly=FALSE,
                           nonSplitOnly=FALSE,
                           primaryOnly=FALSE,
                           ignoreDup=FALSE,
                           
                           # strandness
                           strandSpecific=0,
                           
                           # exon-exon junctions
                           juncCounts=TRUE,
                           genome=NULL,
                           
                           # parameters specific to paired end reads
                           isPairedEnd=TRUE,
                           requireBothEndsMapped=FALSE,
                           checkFragLength=TRUE,
                           minFragLength=10,
                           maxFragLength=10000,
                           countChimericFragments=FALSE,	
                           autosort=TRUE,
                           
                           # miscellaneous
                           nthreads=1,
                           maxMOp=10,
                           reportReads=FALSE)


write.table(x=data.frame(fc_human$annotation, stringsAsFactors=FALSE),file="featureCounts_output/fc_run_all_humanANNOTAION.txt")
write.table(x=data.frame(fc_human$counts, stringsAsFactors=FALSE),file="featureCounts_output/fc_run_all_humanCOUNTS.txt")
write.table(x=data.frame(fc_human$targets, stringsAsFactors=FALSE),file="featureCounts_output/fc_run_all_humanTARGETS.txt")
write.table(x=data.frame(fc_human$stat, stringsAsFactors=FALSE),file="featureCounts_output/fc_run_all_humanSTAT.txt")
