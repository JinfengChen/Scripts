perl runRepeatMasker160.pl ../input/all.con ../input/OS_MULE_TIRs.fa > log 2> log2 &
perl runRepeatMasker160.pl ../input/Gramene.chr.fa ../input/OBa_MULE_TIRs.fa > log 2> log2 &
perl PackMULE.pl --out Gramene.chr.fa.RepeatMasker.out --fasta ../input/Gramene.chr.fa > log 2> log2 &
perl PackMULE_emboss.pl --out wublast/Gramene.chr.fa.RepeatMasker.out --fasta ../input/Gramene.chr.fa > OBa.log 2> OBa.log2 &
perl PackMULE_align.pl --out wublast/Gramene.chr.fa.RepeatMasker.out --fasta ../input/Gramene.chr.fa > OBa1.log 2> OBa1.log2 &
perl PackMULE_align.pl --out wublast/Gramene.chr.fa.RepeatMasker.out --fasta ../input/Gramene.chr.fa > OBa2.log 2> OBa2.log2 &
perl PackMULE_align.pl --out crossmatch/Gramene.chr.fa.RepeatMasker.out --fasta ../input/Gramene.chr.fa > OBa3.log 2> OBa3.log2 &

perl PackMULE_align.pl --out all.con.RepeatMasker.out --fasta ../input/all.con > OS1.log 2> OS1.log2 &
perl PackMULE_align.pl --out wublast/all.con.RepeatMasker.out --fasta ../input/all.con > OS2.log 2> OS2.log2 &


