for infile in ./*.sam
do
~/bin/samtools view -bS ${infile} > ${infile}.bam
done
