#!bin/bash

for i in {1..1000}
do
    ./poa -read_fasta ../../tools/abPOA/evaluation/data/poa-$i -clustal clustal-$i.out -hb blosum80.mat 2> /dev/null
done
