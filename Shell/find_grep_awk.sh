find ../input/fastq/012/ -name "GN*.Maq.p1.map.pileup.chr05" -exec grep "chromosome05" {} \; | awk '$4 > 0' | less -S

