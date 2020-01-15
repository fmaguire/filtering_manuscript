# just do card canonical
set -e 
#ln -s ../../../data/test_metagenome/test_metagenome.fq .
#diamond makedb --in ../../dbs/card_canonical_aa.faa --db card_diamond

# can use max_target_seqs 1 for diamond despite
# the problem in blastx https://github.com/bbuchfink/diamond/issues/232
#/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" diamond blastx -p 2 --max-target-seqs 1 --db card_diamond --outfmt 6 --out default.top_hits -q test_metagenome.fq
#/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" diamond blastx -p 2 --max-target-seqs 1 --db card_diamond --sensitive --outfmt 6 --out sensitive.top_hits -q test_metagenome.fq
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" diamond blastx -p 2 --max-target-seqs 1 --db card_diamond --more-sensitive --outfmt 6 --out more_sensitive.top_hits -q test_metagenome.fq
