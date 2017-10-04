for f in `ls *.fastq.gz | sed 's/R[12]_001.fastq.gz//g' | sort -u`
do
    ~/bin/hisat2 --un-gz ${f}un.fastq --un-conc-gz ${f}unconc.fastq -p 8 --dta --qc-filter --met-file ${f}met_file.txt -x ../indexes/grch38_tran/genome_tran -1 ${f}R1_001.fastq.gz -2 ${f}R2_001.fastq.gz -S ${f}.sam
done
