#!/bin/bash
set -e

# build the random eskape test metagenome 
#cd random_eskape
#./amrtime generate_training -r 250 -m 150 -c 25 -e MSv3 random_eskape.tsv
#rm output.sam
#cd ..

# build the card and card prevalence protein (homology variant efflux) but NOT 
# rRNA spike in metagenome
#cd card_spike_in
#cat ../../CARD_canonical/nucleotide_fasta_protein_* ../../CARD_prevalence/nucleotide_fasta_protein_* > card_nucleotides.fna
## to replace spaces in card prev names
#sed -i 's/ /_/g' card_nucleotides.fna
#python build_spikin.py
#cd ..

# combine two simulated metagenomes
#cat random_eskape/output.fq card_spike_in/card_spike_in_metagenome.fq > test_metagenome.fq 
#zcat random_eskape/output_labels.tsv > test_labels.tsv
#cat card_spike_in/card_spike_in_labels.tsv >> test_labels.tsv
#cat test_labels.tsv | grep -v -P "na\tna\tna\tna" > test_clean_labels.tsv 

# shuffle them to make it nice and random
# proved to use too much memory so fuck it no shuffling
###echo "shuffling metagenome"
###cat test_metagenome.fq | paste - - - - > test.fq
###shuf --random-source=42 test.fq | tr "\t" "\n" > test_metagenome.fq
###rm test.fq 
###
###echo "shuffling labels"
###shuf --random-source=42 test_labels.tsv > test_labels_shuffled.tsv
###mv test_labels_shuffled.tsv test_labels.tsv

# remove non-AGCT or N characters and replace with N
#perl tidy.pl test_metagenome.fq > test_metagenome_clean.fq
