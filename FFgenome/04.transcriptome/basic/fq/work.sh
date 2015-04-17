perl ../../../../tools/bin/fq_all2std.pl fq2fa Root_R8_L1_1.fq > Root_R8_L1_1.fa
perl ../../../../tools/bin/fq_all2std.pl fq2fa Root_R8_L1_2.fq > Root_R8_L1_2.fa
perl ../../../../tools/bin/fq_all2std.pl fq2fa shoot_2_L1_1.fq > shoot_2_L1_1.fa
perl ../../../../tools/bin/fq_all2std.pl fq2fa shoot_2_L1_2.fq > shoot_2_L1_2.fa
cat Root_R8_L1_1.fa Root_R8_L1_2.fa shoot_2_L1_1.fa shoot_2_L1_2.fa > OBa.RNAseq
