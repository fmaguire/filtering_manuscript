# just do card canonical
set -e 
ln -s ../../../data/test_metagenome/test_metagenome_clean.fq test_metagenome.fq
ln -s ../../dbs/card_canonical_nts.fna .

bwa index card_canonical_nts.fna
#
## default minscore 30
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bwa mem -t 2 card_canonical_nts.fna test_metagenome.fq > default.sam
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bwa mem -t 2 -T 15 card_canonical_nts.fna test_metagenome.fq > low_thresh.sam
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bwa mem -t 2 -T 45 card_canonical_nts.fna test_metagenome.fq > high_thresh.sam
#
for bam in *.bam;
do 
    echo $out
    out=$(echo $bam | sed 's/\.bam/_clean\.bam/')
    samtools view -F 0x04 -b $bam > temp.bam
    samtools sort -n -t NM temp.bam > $out
    rm temp.bam
done

