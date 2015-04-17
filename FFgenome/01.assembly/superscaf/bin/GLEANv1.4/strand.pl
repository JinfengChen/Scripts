#OBR_GLEAN_10013909      2971    chr04   2       Ob04g23830.1,2983,624   Ob04g23840.1,1900,1900

my $gff1=parseGFFpos($ARGV[0]);
my $gff2=parseGFFpos($ARGV[1]);
open IN, "$ARGV[2]" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[3] > 0){
       print "$unit[0]\t$unit[1]\t$unit[2]\t$gff1->{$unit[0]}->[2]\t$unit[3]";
       for(my $i=4;$i<@unit;$i++){
          my @array=split(",",$unit[$i]);
          print "\t$array[0],$gff2->{$array[0]}->[2],$array[1],$array[2]";
       }
       print "\n";
    }else{
       print "$unit[0]\t$unit[1]\t$unit[2]\t$gff1->{$unit[0]}->[2]\t$unit[3]\n";
       
    } 

}
close IN;

#####
sub parseGFFpos
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
            $id=$1 if ($id=~/LOC_(.*)/);
        }
        $hash{$id}=[$unit[3],$unit[4],$unit[6]]; 
    }
}
close IN;
return \%hash;
}

