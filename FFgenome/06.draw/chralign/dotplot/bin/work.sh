perl sumBlockV2.pl --block chr01.blocks
perl checkmissinblock.pl --bed ../input/Os.bed.chr/Os01.bed --blockgene chr01.blocks.geneinblock --blast chr01.blast > chr01.check

echo "run chromosome one by one"
perl runChrAll.pl > log 2> log2 &

perl runChrAll.pl --outdir chrmcscan_final8_3_last > log 2> log2 &

perl runChrOryza.pl --outdir chrmcscan_oryza > log 2> log2 &

perl runChrAll.pl --outdir chrmcscan_OBav3 > log 2> log2 &

perl runChrAll.pl --outdir chrmcscan_OBav4 > log 2> log2 &

perl runChrtwo.pl --outdir 2wayv1.4 > log 2> log2 &

perl runChrtwo.pl --outdir 2wayOg > log 2> log2 &

perl runChrtwo.pl --outdir 2wayOi > log 2> log2 &

perl runMcscan.pl --project grass > log 2> log2 &


