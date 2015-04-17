echo "get block inf"
perl GetChrBlockOSREF.pl --align ../input/2wayOi --qrygff ../input/Oi.bed --refgff ../input/Os.bed > TigrVSOi.inf
perl GetChrBlock.pl --align ../input/2wayv1.4 --qrygff ../input/Os.bed --refgff ../input/Ob.bed > TigrVSGramenev1.4.inf
perl GetChrBlock.pl --align ../input/2wayOg --qrygff ../input/Os.bed --refgff ../input/Og.bed > TigrVSOg.inf

echo "get data"
perl GetRegionPipe4ACT.pl --block TigrVSGramenev1.4.inf > log 2> log2 &
perl GetRegionPipe4ACTOSREF.pl --block TigrVSOi.inf > log 2> log2 &

echo "draw graph"
perl drawblockpipe.pl --block graph.inf > log 2> log2 
perl drawblockpipe.pl --block chr04.inf > log 2> log2
perl drawblockpipe.pl --block heter.inf > log 2> log2 &


