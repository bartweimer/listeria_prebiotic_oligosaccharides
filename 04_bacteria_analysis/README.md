# Bacteria Samples

This folder contains codes relevant to the analysis of the listeria-infected human cells. All reads that did not align to the human genome were output in to a fastq file. These reads were then aligned to the listeria genome (GCA_000196035.1_ASM19603v1) using bowtie2. Differential expression was performed on all samples that were infected with listeria (i.e., control cells with no listeria infection were excluded, because approximately no reads aligned to the listeria genome). 
