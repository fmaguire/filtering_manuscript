#/bin/bash
# because we want 100 nt mininmum (we use 150 as our minimum coverage in the 
/usr/bin/time -o time.tsv -a -f "%C\t%e\t%M" orfm -m 96 ../../../data/test_metagenome/test_metagenome.fq > test_metagenome_aa.faa
