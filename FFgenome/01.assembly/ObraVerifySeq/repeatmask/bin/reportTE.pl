## read in a gff3 file of TE annotation, and scaffold file containing scaffold length

## we report name length and TE% in the result file. 

die "Usage: perl reportTE.pl .out.gff lengthfile > log &" if (@ARGV < 1);

my %scaflen;
open IN, "$ARGV[1]" or die "$!";
       while (<IN>){
           chomp $_;
           my @unit=split("\t",$_);
           $scaflen{$unit[0]}=$unit[1];      
       }
close IN;
my %scafte;
open IN, "$ARGV[0]" or die "$!";
       while (<IN>){
           chomp $_;
           my @unit=split("\t",$_);
           #print "$unit[3]\t$unit[4]\n";
           if (exists $scafte{$unit[0]}){
                my $refarray=$scafte{$unit[0]};
                push (@$refarray,"$unit[3]\t$unit[4]");
                $scafte{$unit[0]}=$refarray;
           }else{
                my @pos;
                push (@pos, "$unit[3]\t$unit[4]");
                $scafte{$unit[0]}=\@pos;
           }
       }
close IN;
open OUT1, ">TEcontent" or die "can not open my outfile1";
foreach (keys %scaflen){
      my $scaf=$_;
      if (exists $scafte{$scaf}){
             my $start=1;
             my $end=$scaflen{$scaf}+1;
             my $tepercent=sumte($start,$end,$scafte{$scaf}); 
             print OUT1 "$scaf\t$scaflen{$scaf}\t$tepercent\n";
      }else{
             print OUT1 "$scaf\t$scaflen{$scaf}\t0\n";
      }
}
close OUT1;


sub sumte{
my ($start,$end,$scafte)=@_;
my $scaflen=$end-$start;
my $telen;
foreach (@$scafte){
     my @unit=split("\t",$_);
     if ($unit[0] > $start and $unit[1] < $end){
          $telen+=$unit[1]-$unit[0];
     }elsif($unit[0] < $start and $unit[1] > $start){
          $telen+=$unit[1]-$start;
     }elsif($unit[0] < $end and $unit[1] > $end){
          $telen+=$end-$unit[1];
     }
}
my $percent=$telen/$scaflen;
return $percent;
}
