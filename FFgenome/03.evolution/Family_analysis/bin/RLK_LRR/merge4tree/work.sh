cp ../rice_RLK_LRR/*.LRR*.fa ./
cp ../OBa/type/OBa.LRR*.fa ./
cp ../OBa/type/OBa.Blast.LRR*.fa ./
cp ../SB/type/SB.Blast.LRR*.fa ./
./merge.sh
rm Rice.LRR* SB.Blast.LRR* OBa.Blast.LRR*
perl Drawdir.pl > draw.log 2> draw.log2 &

