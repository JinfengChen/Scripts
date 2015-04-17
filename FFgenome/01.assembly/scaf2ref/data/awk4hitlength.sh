awk '{len=$8-$7;print len}' scaf2chr04blastm8 > awktest
wc -l awktest
