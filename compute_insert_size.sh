for i in 1 2 3 4 5 6 7 8 9 10 
do samtools view $i/accepted_hits.bam|head -n 200 | awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print "n="N", mean="M", stdev="sqrt ((S2-M*M*N)/(N-1))}' 

done
