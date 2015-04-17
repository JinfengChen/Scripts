perl CompareRPKM.pl --bed ../input/rice.mRNA.bed --rpkm1 ../input/RAP3.shoot.rpkm.tophat --rpkm2 ../input/RAP3.shoot.rpkm.tophat --project Rice


perl CompareRPKMorth.pl --ortholog ../input/OB-OS.orth --rpkm1 ../input/FF.shoot.rpkm.tophat --rpkm2 ../input/RAP3.shoot.rpkm.tophat --project FF_Rice

perl Expression2TE.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.repeat.chr.gff --rpkm ../input/OBa.all.tophat.gff.expr --windows 1000 --project OBa_1000 > log 2> log2 &

perl TEdistr2species.pl --ortholog ../input/OB-OS.orth --bed1 ../input/OBa.mRNA.bed --gff1 ../input/OBa.repeat.chr.gff --bed2 ../input/rice.mRNA.bed --gff2 ../input/IRGSP.build5.RepeatMasker.out.gff --windows 500 --project FF_rice_500 &

perl splitRepeatByMethy.pl --methyCG ../input/rice.REPEAT.CG.level.part --methyCHG ../input/rice.REPEAT.CHG.level.part --methyCHH ../input/rice.REPEAT.CHH.level.part -gff ../input/IRGSP.build5.RepeatMasker.out.gff --project rice > log 2> log2 &
perl splitRepeatByMethy.pl --methyCG ../input/OBa.REPEAT.CG.level.part --methyCHG ../input/OBa.REPEAT.CHG.level.part --methyCHH ../input/OBa.REPEAT.CHH.level.part -gff ../input/OBa.all.manual.TE.chr.gff --project OBa > log 2> log2 &
perl splitRepeatByMethy.pl --methyCG ../input/FF.REPEAT.CG.level.part --methyCHG ../input/FF.REPEAT.CHG.level.part --methyCHH ../input/FF.REPEAT.CHH.level.part -gff ../input/FF.repeat.chr.gff --project FF > log 2> log2 &


perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ../input/RAP3.shoot.rpkm.tophat --project rice > log 2> log2 &
perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.repeat.chr.gff --rpkm ../input/OBa.all.tophat.gff.expr --project OBa > log 2> log2 &

perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ../input/RAP3.shoot.rpkm.tophat --me ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status --project rice > log 2> log2 &
perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ../input/OBa.all.tophat.gff.expr --me ../input/OBa.all.manual.TE.chr.gff.Me.status --project OBa > log 2> log2 &
perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/FF.repeat.chr.gff --rpkm ../input/OBa.all.tophat.gff.expr --me ../input/FF.repeat.chr.gff.Me.status --project FF > log 2> log2 &

perl orthTE2Expression.pl --orth ../input/OB-OS.orth --inf1 ./Distance2Expression/0.1/1k/OBa_Me.gene.closestTE.inf --inf2 ./Distance2Expression/0.1/1k/rice_Me.gene.closestTE.inf --project orth_OBa_rice > log 2> log2 &

perl Expression2TE.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ../input/OBa.all.tophat.gff.expr --windows 50000 --project OBa_50000 > log 2> log2 &
perl Expression2TE.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ../input/RAP3.shoot.rpkm.tophat --windows 50000 --project rice_50000 > log 2> log2 &


perl summaryGFF4TE.pl --gff ../input/OBa.all.manual.TE.chr.gff --bed ../input/OBa.mRNA.bed --status ../input/OBa.all.manual.TE.chr.gff.Me.status > log 2> log2 &
perl summaryGFF4TE.pl --gff ../input/IRGSP.build5.RepeatMasker.out.gff --bed ../input/rice.mRNA.bed  --status ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status > log 2> log2 &


perl RPKMwindowbar.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_RPKM_bar --gff ../input/RAP3.gff3.nr.gff.chr --rpkm ../input/RAP3.shoot.rpkm.tophat > log 2> log2 &
perl RPKMwindowbar.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_RPKM_bar --gff ../input/OBa.all.gff.chr --rpkm ../input/OBa.all.tophat.gff.expr > log 2> log2 &

perl splitTEfile.pl --gff ../input/OBa.all.manual.TE.chr.gff
perl splitTEfile.pl --gff ../input/IRGSP.build5.RepeatMasker.out.gff

perl RPKMnormalization.pl --rpkm1 ../input/OBa.all.tophat.gff.expr --rpkm2 ../input/RAP3.shoot.rpkm.tophat --project normalization > log 2> log2 &

perl CompareRPKMorth.pl --ortholog ../input/OB-OS.orth --rpkm1 normalization.1.4r --rpkm2 normalization.2.4r --project OBa_rice > log 2> log2 &


perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm normalization.1.4r --me ../input/OBa.all.manual.TE.chr.gff.Me.status --project OBa > log 2> log2 &
perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm normalization.2.4r --me ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status --project rice > log 2> log2 &

perl orthTE2Expression.pl --orth ../input/OB-OS.orth --inf1 OBa.gene.closestTE.inf --inf2 rice.gene.closestTE.inf --project orth_OBa_rice > log 2> log2 &

perl Expression2TE.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ./normalization/normalization.1.4r --windows 50000 --project OBa_50000 > log 2> log2 &
perl Expression2TE.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm normalization/normalization.2.4r --windows 50000 --project rice_50000 >log 2>log2 &

perl CompareRPKMorth.pl --ortholog ../input/OB-OS.orth --rpkm1 ./normalization/normalization.1.4r --rpkm2 ./normalization/normalization.2.4r --project OBa_rice > log 2> log2 &


perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ./normalization/normalization.1.4r --me ../input/OBa.all.manual.TE.chr.gff.Me.status --type all --project OBa_all >log 2>log2 &
perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ./normalization/normalization.1.4r --me ../input/OBa.all.manual.TE.chr.gff.Me.status --type up --project OBa_5UTR > log 2> log2 &
perl Distance2RPKM.pl --bed ../input/OBa.mRNA.bed --gff ../input/OBa.all.manual.TE.chr.gff --rpkm ./normalization/normalization.1.4r --me ../input/OBa.all.manual.TE.chr.gff.Me.status --type down --project OBa_3UTR > log 2> log2 &

perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ./normalization/normalization.2.4r --me ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status --type all --project rice_all > rice_all.log 2> rice_all.log2 &
perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ./normalization/normalization.2.4r --me ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status --type up --project rice_5UTR > rice_5UTR.log 2> rice_5UTR.log2 &
perl Distance2RPKM.pl --bed ../input/rice.mRNA.bed --gff ../input/IRGSP.build5.RepeatMasker.out.gff --rpkm ./normalization/normalization.2.4r --me ../input/IRGSP.build5.RepeatMasker.out.gff.Me.status --type down --project rice_3UTR > rice_3UTR.log 2> rice_3UTR.log2 &

perl orthTE2Expression.pl --orth ../input/OB-OS.orth --inf1 ./Distance2Expression/allTE/0.1/normalized/type/1kb/OBa_5UTR.gene.closestTE.inf --inf2 ./Distance2Expression/allTE/0.1/normalized/type/1kb/rice_5UTR.gene.closestTE.inf --project orth_OBa_rice_5UTR >log 2>log2 &
perl orthTE2Expression.pl --orth ../input/OB-OS.orth --inf1 ./Distance2Expression/allTE/0.1/normalized/type/1kb/OBa_3UTR.gene.closestTE.inf --inf2 ./Distance2Expression/allTE/0.1/normalized/type/1kb/rice_3UTR.gene.closestTE.inf --project orth_OBa_rice_3UTR > log 2> log2 &
perl orthTE2Expression.pl --orth ../input/OB-OS.orth --inf1 ./Distance2Expression/allTE/0.1/normalized/type/1kb/OBa_all.gene.closestTE.inf --inf2 ./Distance2Expression/allTE/0.1/normalized/type/1kb/rice_all.gene.closestTE.inf --project orth_OBa_rice_all > log 2> log2 &


