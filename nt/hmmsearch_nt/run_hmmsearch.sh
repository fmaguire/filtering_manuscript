#!/bin/bash
set -e 

# dump to family_fasta
# singles get converted to stockholm notation
# multiple member families get aligned with mafft --auto and converted to sto
# hmmbuild with default settings
# cat *.hmm
# hmmpress

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" hmmsearch --cpu 2 --noali --notextw --tblout nt_hmmsearch_default.tbl nt_card.hmm test_metagenome.fna 

tail -n +4 nt_hmmsearch_default.tbl | sort -T ~/oveflow_data/tmp_sort -r -k 1,5 > nt_hmmsearch_default_sorted.tbl
python sort_hmmsearch_by_evalue.py 
