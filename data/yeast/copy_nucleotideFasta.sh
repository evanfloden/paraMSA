#!/bin/bash
gene_list="/users/cn/jchang/data/yeast/2008_science-yeast1502/geneList_YeastToL.txt" #853 sets

nucleotide_fasta_p="../nucleotide_fasta"
cmd="t_coffee -other_pg seq_reformat"

cat $gene_list|while read set
do
  echo "process $set"
  $cmd -in $nucleotide_fasta_p/prank/$set.fasta -action +translate -output fasta_aln|sed s/o/-/g > prank/$set.fasta
  $cmd -in $nucleotide_fasta_p/sate/$set.fasta -action +translate -output fasta_aln|sed s/o/-/g > sate/$set.fasta
done
