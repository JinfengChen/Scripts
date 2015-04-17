bowtie-build ../input/LargeVerifySeq.fas LargeVerifySeq >log 2>log2
tophat -r 50 LargeVerifySeq ../fq/Root_R8_L1_1.fq ../fq/Root_R8_L1_2.fq > log 2> log2 &

bowtie-build ../input/OBa.all.fa OBa.all > log 2> log2 &
perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -1 ../fq/Root_1 -2 ../fq/Root_2 -p OBaTest

perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -1 ../fq/Root_R8_L1_1.fq -2 ../fq/Root_R8_L1_2.fq -p OBaRoot > log 2> log2 &

