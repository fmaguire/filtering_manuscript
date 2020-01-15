#!/bin/bash
set -e

#ln -s ../../../data/test_metagenome/test_metagenome.fq .
#
#ln -s ../../dbs/card_canonical_aa.faa .
#paladin index -r3 card_canonical_aa.faa

#/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" paladin align -f 100 -t 2 -T 5 -C card_canonical_aa.faa test_metagenome.fq | samtools view -Sb > min5.bam
#/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" paladin align -f 100 -t 2 -T 15 -C card_canonical_aa.faa test_metagenome.fq | samtools view -Sb > min15.bam
#/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" paladin align -f 100 -t 2 -T 25 -C card_canonical_aa.faa test_metagenome.fq | samtools view -Sb > min25.bam
#

for bam in *.bam;
do 
    out=$(echo $bam | sed 's/\.bam/_clean\.bam/')
    echo $out
    samtools view -F 0x04 -b $bam > temp.bam
    samtools sort -n -t NM temp.bam > $out
    rm temp.bam
done

