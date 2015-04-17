cat tigr.all.cds.final.clean.fa Sb.final.clean.fa > OsSb.cds.fa
echo "blast bradi cds to Os and Ob cds using blastn"
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -i Bradi_1.0.cds.final.clean.fa -d OsSb.cds.fa > log 2> log2 &


