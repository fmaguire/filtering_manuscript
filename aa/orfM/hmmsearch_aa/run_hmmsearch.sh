#!/bin/bash
set -e 
#fastq_to_fasta -i ../../../data/test_metagenome/test_metagenome_clean.fq -o ../../../data/test_metagenome/test_metagenome.fna 
ln -sfn ../test_metagenome_aa.faa .
#bash prepare_hmms.sh

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" hmmsearch --cpu 2 --noali --notextw --tblout aa_hmmsearch_default.tbl aa_card.hmm test_metagenome_aa.faa

tail -n +4 aa_hmmsearch_default.tbl | sort -T ~/overflow_data/tmp_sort -r -k 1,5 > aa_hmmsearch_default_sorted.tbl
