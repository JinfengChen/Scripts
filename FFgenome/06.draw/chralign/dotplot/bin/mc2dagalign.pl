##Convert the result of mcscan to DAGalign, So I can use perl script to draw dotplot 
##use .gff and .align output 4dotplot
##Write by Chenjinfeng on 20090929
##Useage: perl mc2dagalign.pl os_sb ## os_sb is the prefix for aligns and gff file, for example os_sb.gff and os_sb.aligns



chomp $ARGV[0];
my $name=$ARGV[0];
#print "$name\n";

######read in gff file and store infor in hash## 
my %chromosome;
my %position;
open GFF, "$name.gff" or die "can not open my gff file";
while(<GFF>){
    chomp $_;
    my ($chr,$gene,$five,$three)= split("\t",$_);
    if($chr=~/(\d+)/){$chr=$1};
    #print "$chr\n";
    $chromosome{$gene}="$chr";
    $position{$gene}="$five\t$three";
}
close GFF;

#######################################3
my $outfile="$name"."4dotplot";
open OUT, ">$outfile" or die "can not open my out file"; 
open ALIGN, "$name.aligns" or die "can not open my align file";
my $i=1;
while($i<=11){
$i++;
<ALIGN>;
}
while(<ALIGN>){
      chomp $_;
      unless($_=~/\#\# Alignment/){
          my @unit=split("\t",$_); 
          #print "$unit[1]\t$unit[2]\n";
          my $locus1=$unit[1];
          my $locus2=$unit[2];
          my $evalue=$unit[3];
          print OUT "$chromosome{$locus1}\t$locus1\t$position{$locus1}\t$chromosome{$locus2}\t$locus2\t$position{$locus2}\t$evalue\n";
      }else{
          print OUT "$_\n";
      }
}
close ALIGN;
close OUT;
