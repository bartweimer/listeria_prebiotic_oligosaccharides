setwd("/Users/bcweimer/Desktop/Po_seq/hisat_output/hisat_unconc/listeria/bowtie_bam/")
library(Rsubread)

fls <- dir(".", "bam$")
fls
# This Run:
# [1] "AACAACCA_S7_L005_unconc_bac.sam.bam"  "AACCGAGA_S6_L005_unconc_bac.sam.bam" 
# [3] "AACGCTTA_S5_L005_unconc_bac.sam.bam"  "AAGGTACA_S23_L007_unconc_bac.sam.bam"
# [5] "ACACAGAA_S24_L007_unconc_bac.sam.bam" "ACAGCAGA_S25_L007_unconc_bac.sam.bam"
# [7] "ACCTCCAA_S26_L007_unconc_bac.sam.bam" "ACGCTCGA_S27_L007_unconc_bac.sam.bam"
# [9] "ACGTATCA_S20_L007_unconc_bac.sam.bam" "ACTATGCA_S21_L007_unconc_bac.sam.bam"
# [11] "AGAGTCAA_S22_L007_unconc_bac.sam.bam" "AGATCGCA_S28_L008_unconc_bac.sam.bam"
# [13] "AGCAGGAA_S29_L008_unconc_bac.sam.bam" "AGTCACTA_S30_L008_unconc_bac.sam.bam"
# [15] "ATCCTGTA_S31_L008_unconc_bac.sam.bam" "ATTGAGGA_S32_L008_unconc_bac.sam.bam"
# [17] "CAACCACA_S33_L008_unconc_bac.sam.bam" "CATCAAGT_S8_L005_unconc_bac.sam.bam" 
# [19] "GACTAGTA_S34_L008_unconc_bac.sam.bam"

fc_listeria<-featureCounts(files=fls,
                             # annotation
                             annot.inbuilt=NULL,
                             annot.ext="GCA_000196035.1_ASM19603v1_genomic.gtf",
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


write.table(x=data.frame(fc_listeria$annotation, stringsAsFactors=FALSE),file="fc_listeriaANNOTAION.txt")
write.table(x=data.frame(fc_listeria$counts, stringsAsFactors=FALSE),file="fc_listeriaCOUNTS.txt")
write.table(x=data.frame(fc_listeria$counts_junction, stringsAsFactors=FALSE),file="fc_listeriaCOUNTJUNC.txt")
write.table(x=data.frame(fc_listeria$targets, stringsAsFactors=FALSE),file="fc_listeriaTARGETS.txt")
write.table(x=data.frame(fc_listeria$stat, stringsAsFactors=FALSE),file="fc_listeriaSTAT.txt")

# 
# ==========     _____ _    _ ____  _____  ______          _____  
# =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
# =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
#   ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
#   ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
#   ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
#   Rsubread 1.24.0
# 
# //========================== featureCounts setting ===========================\\
# ||                                                                            ||
#   ||             Input files : 19 BAM files                                     ||
#   ||                           P AACAACCA_S7_L005_unconc_bac.sam.bam            ||
#   ||                           P AACCGAGA_S6_L005_unconc_bac.sam.bam            ||
#   ||                           P AACGCTTA_S5_L005_unconc_bac.sam.bam            ||
#   ||                           P AAGGTACA_S23_L007_unconc_bac.sam.bam           ||
#   ||                           P ACACAGAA_S24_L007_unconc_bac.sam.bam           ||
#   ||                           P ACAGCAGA_S25_L007_unconc_bac.sam.bam           ||
#   ||                           P ACCTCCAA_S26_L007_unconc_bac.sam.bam           ||
#   ||                           P ACGCTCGA_S27_L007_unconc_bac.sam.bam           ||
#   ||                           P ACGTATCA_S20_L007_unconc_bac.sam.bam           ||
#   ||                           P ACTATGCA_S21_L007_unconc_bac.sam.bam           ||
#   ||                           P AGAGTCAA_S22_L007_unconc_bac.sam.bam           ||
#   ||                           P AGATCGCA_S28_L008_unconc_bac.sam.bam           ||
#   ||                           P AGCAGGAA_S29_L008_unconc_bac.sam.bam           ||
#   ||                           P AGTCACTA_S30_L008_unconc_bac.sam.bam           ||
#   ||                           P ATCCTGTA_S31_L008_unconc_bac.sam.bam           ||
#   ||                           P ATTGAGGA_S32_L008_unconc_bac.sam.bam           ||
#   ||                           P CAACCACA_S33_L008_unconc_bac.sam.bam           ||
#   ||                           P CATCAAGT_S8_L005_unconc_bac.sam.bam            ||
#   ||                           P GACTAGTA_S34_L008_unconc_bac.sam.bam           ||
#   ||                                                                            ||
#   ||             Output file : ./.Rsubread_featureCounts_pid36202               ||
#   ||                 Summary : ./.Rsubread_featureCounts_pid36202.summary       ||
#   ||              Annotation : GCA_000196035.1_ASM19603v1_genomic.gtf (GTF)     ||
#   ||       Junction Counting : <output_file>.jcounts                            ||
#   ||      Dir for temp files : .                                                ||
#   ||                                                                            ||
#   ||                 Threads : 1                                                ||
#   ||                   Level : feature level                                    ||
#   ||              Paired-end : yes                                              ||
#   ||         Strand specific : no                                               ||
#   ||      Multimapping reads : counted (as fractions)                           ||
#   || Multi-overlapping reads : counted                                          ||
#   ||       Overlapping bases : 10                                               ||
#   ||       Overlapping bases : 0.0%                                             ||
#   ||                                                                            ||
#   ||          Chimeric reads : not counted                                      ||
#   ||        Both ends mapped : not required                                     ||
#   ||         Fragment length : 10 - 10000                                       ||
#   ||                                                                            ||
#   \\===================== http://subread.sourceforge.net/ ======================//
#     
#     //================================= Running ==================================\\
#     ||                                                                            ||
#       || Load annotation file GCA_000196035.1_ASM19603v1_genomic.gtf ...            ||
#       ||    Features : 1262                                                         ||
#       ||    Meta-features : 1151                                                    ||
#       ||    Chromosomes/contigs : 1                                                 ||
#       ||                                                                            ||
#       || Process BAM file AACAACCA_S7_L005_unconc_bac.sam.bam...                    ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 2858781                                               ||
#       ||    Successfully assigned fragments : 15 (0.0%)                             ||
#       ||    Running time : 0.39 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AACCGAGA_S6_L005_unconc_bac.sam.bam...                    ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 6796587                                               ||
#       ||    Successfully assigned fragments : 11 (0.0%)                             ||
#       ||    Running time : 0.88 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AACGCTTA_S5_L005_unconc_bac.sam.bam...                    ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 8249696                                               ||
#       ||    Successfully assigned fragments : 5 (0.0%)                              ||
#       ||    Running time : 0.96 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AAGGTACA_S23_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 31583393                                              ||
#       ||    Successfully assigned fragments : 1147413 (3.6%)                        ||
#       ||    Running time : 3.80 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACACAGAA_S24_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 5040834                                               ||
#       ||    Successfully assigned fragments : 100131 (2.0%)                         ||
#       ||    Running time : 0.44 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACAGCAGA_S25_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 22136322                                              ||
#       ||    Successfully assigned fragments : 460478 (2.1%)                         ||
#       ||    Running time : 3.01 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACCTCCAA_S26_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 4205561                                               ||
#       ||    Successfully assigned fragments : 67727 (1.6%)                          ||
#       ||    Running time : 0.52 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACGCTCGA_S27_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 5452099                                               ||
#       ||    Successfully assigned fragments : 44917 (0.8%)                          ||
#       ||    Running time : 0.86 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACGTATCA_S20_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 26792666                                              ||
#       ||    Successfully assigned fragments : 514464 (1.9%)                         ||
#       ||    Running time : 4.10 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ACTATGCA_S21_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 34531603                                              ||
#       ||    Successfully assigned fragments : 468756 (1.4%)                         ||
#       ||    Running time : 4.47 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AGAGTCAA_S22_L007_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 11367644                                              ||
#       ||    Successfully assigned fragments : 207531 (1.8%)                         ||
#       ||    Running time : 1.76 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AGATCGCA_S28_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 15728023                                              ||
#       ||    Successfully assigned fragments : 113078 (0.7%)                         ||
#       ||    Running time : 2.31 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AGCAGGAA_S29_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 22299059                                              ||
#       ||    Successfully assigned fragments : 400436 (1.8%)                         ||
#       ||    Running time : 3.53 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file AGTCACTA_S30_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 17699889                                              ||
#       ||    Successfully assigned fragments : 538365 (3.0%)                         ||
#       ||    Running time : 2.79 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ATCCTGTA_S31_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 8210726                                               ||
#       ||    Successfully assigned fragments : 196562 (2.4%)                         ||
#       ||    Running time : 1.35 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file ATTGAGGA_S32_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 10648965                                              ||
#       ||    Successfully assigned fragments : 579 (0.0%)                            ||
#       ||    Running time : 1.62 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file CAACCACA_S33_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 14009392                                              ||
#       ||    Successfully assigned fragments : 568 (0.0%)                            ||
#       ||    Running time : 1.11 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file CATCAAGT_S8_L005_unconc_bac.sam.bam...                    ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 27724311                                              ||
#       ||    Successfully assigned fragments : 50 (0.0%)                             ||
#       ||    Running time : 2.02 minutes                                             ||
#       ||                                                                            ||
#       || Process BAM file GACTAGTA_S34_L008_unconc_bac.sam.bam...                   ||
#       ||    Paired-end reads are included.                                          ||
#       ||    Assign fragments (read pairs) to features...                            ||
#       ||    Total fragments : 6968463                                               ||
#       ||    Successfully assigned fragments : 291 (0.0%)                            ||
#       ||    Running time : 0.48 minutes                                             ||
#       ||                                                                            ||
#       ||                 Found 0 junctions in all the input files.                  ||
#       ||                                                                            ||
#       ||                         Read assignment finished.                          ||
#       ||                                                                            ||
#       \\===================== http://subread.sourceforge.net/ ======================//
