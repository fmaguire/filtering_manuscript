#!/bin/bash
set -e

#mkdir -p card_fasta
#python split_individual.py ../../dbs/card_canonical_nts.fna

#ln -s ../../dbs/card_canonical_nts.fna .
#ln -s ../../../data/test_metagenome/test_metagenome.fq .

mkdir -p filters
mkdir -p output

for k in 5 9 15 25
do
    ./biobloommaker -k $k -p card_${k} card_canonical_nts.fna
    mv card_${k}.* filters

    # if score is 1 the best hit is used and the score appended
    /usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" ./biobloomcategorizer -t 2 -i -f filters/card_${k}.bf -s 1 --fa -g -p output/${k} test_metagenome.fq 


    rm output/${k}_noMatch.fa.gz output/${k}_multiMatch.fa.gz
done

