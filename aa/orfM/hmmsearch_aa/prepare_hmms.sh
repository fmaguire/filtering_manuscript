#!/bin/bash
set -e 

# dump to family_fasta
# singles get converted to stockholm notation
# multiple member families get aligned with mafft --auto and converted to sto
# hmmbuild with default settings
# cat *.hmm
# hmmpress
# dump fasta

source ~/miniconda3/etc/profile.d/conda.sh
conda activate AMRtime_analysis

ln -sfn ../../../scripts/utils.py .
python dump_to_gene_family_fasta.py
#
#mkdir -p family_fasta_aa/hmms/{single,group}
#
#cd family_fasta_aa
#
#for fasta_file in $(find -type f)
#do
#    seqs=$(grep -c "^>" $fasta_file)
#    
#    if [ "$seqs" -gt 1 ]
#    then
#        cp $fasta_file hmms/group
#    fi
#
#    if [ "$seqs" -eq 1 ]
#    then
#        cp $fasta_file hmms/single
#    fi
#done
#
#cd hmms/single
#for seq in $(find -type f)
#do
#    esl-reformat stockholm $seq > ${seq}.sto
#    hmmbuild ${seq}.hmm ${seq}.sto
#done
#
#cd ../group
#for seq in $(find -type f)
#do
#    mafft --auto $seq > ${seq}.afa
#    esl-reformat stockholm ${seq}.afa > ${seq}.sto
#    hmmbuild ${seq}.hmm ${seq}.sto
#done

#cd family_fasta_aa/hmms/group
#cd ..
#find -type f -name '*.hmm' -exec cat {} \; > ../../aa_card.hmm
#cd ../..
#hmmpress aa_card.hmm

