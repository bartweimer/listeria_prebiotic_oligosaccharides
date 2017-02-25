for f in `ls *.fastq.gz | sed 's/.[12].fastq.gz//g' | sort -u`
	do
	~/bin/bowtie2 --met-file ${f}_bac_met.txt --un-conc-gz ${f}_bac_unconc%.fastq.gz -x index1/listeria -1 ${f}.1.fastq.gz -2 ${f}.2.fastq.gz -S ${f}_bac.sam
done
