###get gc content for gcdraw using sliding windos and get chromosome infor
###write by chenjinfeng,modify on 20090928

$/=">";
#open IN, "EC611973.txt";
open IN, "/home/seqlib/sorghum/Sorbi1_assembly_scaffolds.fasta" or die "can not open my infile";
open NU, ">sorghumChrStatistic.txt" or die "can not open my rice static";
while (<IN>){
     my @unit=split("\n");
     my $head=shift @unit;
     my $seq=join("",@unit);
     $seq=~s/\>//;
     if ($head=~/chr/){
        open OUT, ">$head\.gc" or die "can not open gc out file";
        my $length=length $seq;
        print NU "$head\t$length\n";
        my $step=100000;
        my $win=500000;
        my $run=int (($length-$win)/$step)+1;
        for (my $i=0;$i<=$run;$i++){
            my $start =$i*$step;
            #print "$start\n";
            my $subseq=substr($seq,$start,$win);
           
            my $mid=$start+$win/2;
            my $gc=estimateGC($subseq);
            print OUT "$mid\t$gc\n";
        }
        close OUT;
     }
}
close IN;
close NU;


sub estimateGC{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
#my $t=$seq=~tr/Tt/Tt/g;
my $c=$seq=~tr/Cc/Cc/;
my $g=$seq=~tr/Gg/Gc/;
my $gc=($g+$c)/$length;
return $gc;

}
