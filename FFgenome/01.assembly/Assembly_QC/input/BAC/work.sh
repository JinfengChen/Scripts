echo "check gap in BAC sequence"
perl /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/checkgap.pl --fasta BAC.fasta
echo "repeatmasker BAC sequence with FFTE lib"
perl /home/jfchen/FFproject/tools/bin/runRepeatMasker160.pl BAC.fasta /home/jfchen/FFproject/FFgenome/02.annotation/01.repeat/1.repeat/manual/TElib-FF > log 2> log2 &

