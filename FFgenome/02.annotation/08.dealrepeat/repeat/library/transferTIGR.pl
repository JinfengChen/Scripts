#!/usr/bin/perl

### transfer TIGRrepeat to lib file that needed by repeatmasker
### >Name#class/subclass/family
### Usage: perl transferTIGR.pl TIGRrepeat OUTFILE > log &;

die "Usage: perl transferTIGR.pl TIGRrepeat OUTFILE > log &\n" if (@ARGV < 2);

$/=">";
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[1]" or die "$!";
while(<IN>){
     chomp $_;
     next if (length $_ < 2); 
     my @unit=split("\n",$_);
     my $head=shift @unit;
     my $seq =join("\n",@unit);
     my @array=split(" ",$head);
     my $type=substr($array[0],4,7);
     my $name=pop @array;
     #print "$type\n";
     my $class;
     my $subclass;
     if ($type=~/TERT/){
        if ($type=~/001/){
            $class="LTR";
            $subclass="Copia";
        }elsif($type=~/002/){
            $class="LTR";
            $subclass="Gypsy";
        }elsif($type=~/003/){
            $class="LINE";
            $subclass="LINE";
        }elsif($type=~/004/){
            $subclass="SINE";
            $class="SINE";
        }else{
            $class="Other";
            $subclass="Unclassified Retro";
        }
     }elsif($type=~/TETN/){
        $class="DNA";
        if ($type=~/001/){
             $subclass="Ac/Ds";
        }elsif($type=~/002/){
             $subclass="En-Spm";
        }elsif($type=~/003/){
             $subclass="Mutator(MULE)";
        }else{
             $subclass="Unclassified";
        }
     }elsif($type=~/TEMT/){
        $class="DNA";
        if ($type=~/001/){
             $subclass="Tourist";
        }elsif($type=~/002/){
             $subclass="Stowaway";
        }else{
             $subclass="MITE";
        }
     }elsif($type=~/CMCM/){
        $class="Other";
        if($type=~/001/){
             $subclass="Centromere-specific retrotransposon";
        }elsif($type=~/002/){
             $subclass="Centromeric satellite";
        }else{
             $subclass="Unclassified Centromere"; 
        }
     }elsif($type=~/TRTM/){
        $class="Other";
        $subclass="Telomere-associated";
     }elsif($type=~/RGRR/){
        $class="Other";
        $subclass="rDNA";
     }elsif($type=~/OTOT/){
        $class="Other"; 
        $subclass="Unclassified";
     }else{
        $class="Other";        
        $subclass="Unclassified";
     }
     print OUT ">$name#$class/$subclass $array[0]\n$seq\n";
  
}
close OUT;
close IN;
