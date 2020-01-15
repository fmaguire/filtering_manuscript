# just do card canonical
set -e 
ln -sfn ../../../data/test_metagenome/test_metagenome_clean.fq test_metagenome.fq
#
mkdir -p card_db
cd card_db
#
ln -sfn ../../../dbs/card_canonical_nts.fna .

vsearch --cluster_size card_canonical_nts.fna --id 0.90 --msaout MSA.temp

awk '!a[$0]++ {of="./cluster-" ++fc ".msa"; print $0 >> of ; close(of)}' RS= ORS="\n\n" MSA.temp 
rm MSA.temp
cd ..


# low default and high k
groot index -p 2 -i card_db -o card_k5_s64_j090  -k 5 -l 250 -s 64 -j 0.90
groot index -p 2 -i card_db -o card_k5_s64_j099  -k 5 -l 250 -s 64 -j 0.99
groot index -p 2 -i card_db -o card_k5_s128_j090 -k 5 -l 250 -s 128 -j 0.90
groot index -p 2 -i card_db -o card_k5_s128_j099 -k 5 -l 250 -s 128 -j 0.99
groot index -p 2 -i card_db -o card_k5_s256_j090 -k 5 -l 250 -s 256 -j 0.90
groot index -p 2 -i card_db -o card_k5_s256_j099 -k 5 -l 250 -s 256 -j 0.99
echo "k5 done"
groot index -p 2 -i card_db -o card_k7_s64_j090  -k 7 -l 250 -s 64 -j 0.90
groot index -p 2 -i card_db -o card_k7_s64_j099  -k 7 -l 250 -s 64 -j 0.99
groot index -p 2 -i card_db -o card_k7_s128_j090 -k 7 -l 250 -s 128 -j 0.90
root index -p 2 -i card_db -o test_card_k7_s128_j099 -k 7 -l 250 -s 128 -j 0.99 # default
root index -p 2 -i card_db -o test_default -k 7 -l 250 -s 128 -j 0.99 # default
groot index -p 2 -i card_db -o card_k7_s256_j090 -k 7 -l 250 -s 256 -j 0.90
groot index -p 2 -i card_db -o card_k7_s256_j099 -k 7 -l 250 -s 256 -j 0.99
echo "k7 done"
groot index -p 2 -i card_db -o card_k9_s64_j090  -k 9 -l 250 -s 64 -j 0.90
groot index -p 2 -i card_db -o card_k9_s64_j099  -k 9 -l 250 -s 64 -j 0.99
groot index -p 2 -i card_db -o card_k9_s128_j090 -k 9 -l 250 -s 128 -j 0.90
groot index -p 2 -i card_db -o card_k9_s128_j099 -k 9 -l 250 -s 128 -j 0.99
groot index -p 2 -i card_db -o card_k9_s256_j090 -k 9 -l 250 -s 256 -j 0.90
groot index -p 2 -i card_db -o card_k9_s256_j099 -k 9 -l 250 -s 256 -j 0.99
echo "k9 done"

mkdir -p align




for db in card_k* test_default test_card_k*
do
    echo $db
    /usr/bin/time -o test_time.tsv -a -f "%C\t%e\t%M" groot align -p 2 -i $db -f test_metagenome.fq  > test_${db}_default.bam
done

echo "sorting"
for bam in test_*.bam;
do 
    out=$(echo $bam | sed 's/\.bam/_clean\.bam/')
    echo $out
    samtools view -F 0x04 -b $bam > temp.bam
    samtools sort -n -t NM temp.bam > $out
    rm temp.bam
done
