# just do card canonical
set -e 
ln -s ../../../data/test_metagenome/test_metagenome_clean.fq test_metagenome.fq
ln -s ../../dbs/card_canonical_nts.fna .

bowtie2-build card_canonical_nts.fna card_canonical_nts

 default (end-to-end --sensitive
/usr/bin/time -o time2.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --very-fast -x card_canonical_nts test_metagenome.fq > e2e_very_fast.sam
samtools view -S -b e2e_very_fast.sam > e2e_very_fast.bam
rm e2e_very_fast.sam

/usr/bin/time -o time2.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --fast -x card_canonical_nts test_metagenome.fq > e2e_fast.sam
samtools view -S -b e2e_fast.sam > e2e_fast.bam
rm e2e_fast.sam

/usr/bin/time -o time2.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --sensitive -x card_canonical_nts test_metagenome.fq > default.sam
samtools view -S -b default.sam > default.bam
rm default.sam

/usr/bin/time -o time2.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --very-sensitive -x card_canonical_nts test_metagenome.fq > e2e_very_sensitive.sam
samtools view -S -b e2e_very_sensitive.sam > e2e_very_sensitive.bam
rm e2e_very_sensitive.sam

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --very-fast-local -x card_canonical_nts test_metagenome.fq > l_very_fast.sam
samtools view -S -b l_very_fast.sam > l_very_fast.bam
rm l_very_fast.sam

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --fast-local -x card_canonical_nts test_metagenome.fq > l_fast.sam
samtools view -S -b l_fast.sam > l_fast.bam
rm l_fast.sam

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --sensitive-local -x card_canonical_nts test_metagenome.fq > l_sensitive.sam
samtools view -S -b l_sensitive.sam > l_sensitive.bam
rm l_sensitive.sam

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" bowtie2 -p 2 --very-sensitive-local -x card_canonical_nts test_metagenome.fq > l_very_sensitive.sam
samtools view -S -b l_very_sensitive.sam > l_very_sensitive.bam
rm l_very_sensitive.sam

for bam in *.bam;
do 
    echo $bam
    out=$(echo $bam | sed 's/\.bam/_clean\.bam/')
    echo $out
    samtools view -F 0x04 -b $bam > temp.bam
    samtools sort -n -t NM temp.bam > $out
    rm temp.bam
done
