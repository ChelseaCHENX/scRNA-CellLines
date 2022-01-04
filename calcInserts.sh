cd /bgfs/alee/10Xseq/1904_10X/mkfastq/H7NG7DRXX/AL1
outdir=/bgfs/alee/chelsea/projects/10X/CellLine/data/GEO_Submission

for f in AL1_S1_L001_R2_001 AL1_S1_L002_R2_001
do
    cat $f.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $outdir/$f.length.txt
done