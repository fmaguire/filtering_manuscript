#!/bin/bash
set -e 
ln -s ../test_metagenome_aa.faa .
makeblastdb -dbtype prot -in ../../../dbs/card_canonical_aa.faa -out card_canonical

# default word size is 3
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" blastp -query test_metagenome_aa.faa -db card_canonical -out default.out6 -evalue 1e-3 -outfmt 6 -num_threads 2 
# remove all but top hit (due to 'max_target_seqs' kerfuffle)
awk '! a[$1]++' default.out6 > default.top_hits
rm default.out6

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" blastp -query test_metagenome_aa.faa -db card_canonical -out longword.out6 -evalue 1e-3 -outfmt 6 -num_threads 2 -word_size 4
awk '! a[$1]++' longword.out6  > longword.top_hits
rm longword.out6

/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" blastp -query test_metagenome_aa.faa -db card_canonical -out shortword.out6 -evalue 1e-3 -outfmt 6 -num_threads 2 -word_size 2
awk '! a[$1]++' shortword.out6  > shortword.top_hits
rm shortword.out6

