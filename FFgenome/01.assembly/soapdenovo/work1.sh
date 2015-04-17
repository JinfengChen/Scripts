/home/jfchen/bgitraining/denovo/SOAPdenovo_Release1.0/SOAPdenovo pregraph -s chr04.config -d -o chr04 > logpre
echo pregraph Done
/home/jfchen/bgitraining/denovo/SOAPdenovo_Release1.0/SOAPdenovo contig -g chr04 > logcontig
echo contig Done
/home/jfchen/bgitraining/denovo/SOAPdenovo_Release1.0/SOAPdenovo map -s chr04.config  -g chr04 > logmap
echo map Done
/home/jfchen/bgitraining/denovo/SOAPdenovo_Release1.0/SOAPdenovo scaff -g chr04 > logscaff
echo scaff Done
 
