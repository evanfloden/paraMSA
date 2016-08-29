tips=(16)
syms=('0.5')

for tip in ${tips[@]}
do 
  for sym in ${syms[@]}
  do
    rm -rf tips${tip}_${sym}

    mkdir tips${tip}_${sym}

    for i in {001..100}
    do

     #esl-reformat afa /users/cn/jchang/projects/2010-03_ConM-Coffee/results/2013-03-26_bench-simulation/data/tips${tip}/symmetric_${sym}/${i}/seqs/${i}.0400.fa > tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.afa

     #t_coffee -other_pg seq_reformat -in /users/cn/jchang/projects/2010-03_ConM-Coffee/results/2013-03-26_bench-simulation/data/tips${tip}/symmetric_${sym}/${i}/seqs/${i}.0400.fa -output fasta_seq > tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa
   
    cp /users/cn/jchang/projects/2010-03_ConM-Coffee/results/2013-03-26_bench-simulation/data/tips${tip}/asymmetric_${sym}/${i}/seqs/${i}.0400.fa tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa.tmp
    
    t_coffee -other_pg seq_reformat -in tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa.tmp -output fasta_seq -out tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa

    sed -i '/^\s*$/d' tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa

    rm tips${tip}_${sym}/tips${tip}_${sym}_${i}.0400.fa.tmp
    
    done
  done
done

