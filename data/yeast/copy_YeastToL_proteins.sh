#!/bin/bash
gene_list="geneList_YeastToL.txt" #853 sets that have one alignment matching tree of life

path_to_fastas="/users/cn/jchang/projects/2010-03_ConM-Coffee/data/aln/protein/clustalw/"
cmd="t_coffee -other_pg seq_reformat"

cat $gene_list|while read set
do
  echo "process $set"
  $cmd -in $path_to_fastas/$set.fasta -output fasta_seq > protein/$set.fasta
done
