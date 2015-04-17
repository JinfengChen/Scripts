## read in a gff3 file of TE annotation, and scaffold file containing scaffold length
## we report TE% if the length is small than 200kb
## we report TE% and which 200kb windows has a TE% smaller that cut off
## generate two files: One containing scaffold name, length, TE% (TEcontent)
##                     The other containing scaffold name, length, TE% < 10%, start-end on scaffold if > 200k (4FISHprobe)

die "Usage: perl report200kTE.pl windows_size > log &" if (@ARGV < 1);
my $win=$ARGV[0];
my %scaflen;
open IN, "/share/raid12/chenjinfeng/FFgenome/repeat/repeatmasker/data/scaffold" or die "can not open my gff file";
       while (<IN>){
           chomp $_;
           my @unit=split("\t",$_);
           $scaflen{$unit[0]}=$unit[1];      
       }
close IN;
my %scafte;
open IN, "/share/raid12/chenjinfeng/FFgenome/repeat/repeatmasker/output/ass.scafSeq.gapfill3.RepeatMasker.out.gff";
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
open OUT1, ">/share/raid12/chenjinfeng/FFgenome/repeat/repeatmasker/output/TEcontent" or die "can not open my outfile1";
open OUT2, ">/share/raid12/chenjinfeng/FFgenome/repeat/repeatmasker/output/4FISHprobe" or die "can not open my outfile2";
foreach (keys %scaflen){
      my $scaf=$_;
      if (exists $scafte{$scaf}){
         if ($scaflen{$scaf} <= $win){
             my $start=1;
             my $end=$scaflen{$scaf}+1;
             my $tepercent=sumte($start,$end,$scafte{$scaf}); 
             if ($tepercent <= 0.1){
                print OUT1 "$scaf\t$scaflen{$scaf}\t$tepercent\n";  
                #print OUT2 "$scaf\t$scaflen{$scaf}\t$tepercent\n";
             }else{
                print OUT1 "$scaf\t$scaflen{$scaf}\t$tepercent\n";
             }
         }else{
             my $start=1;
             my $end=$scaflen{$scaf}+1;
             my $tepercent=sumte($start,$end,$scafte{$scaf});
             print OUT1 "$scaf\t$scaflen{$scaf}\t$tepercent\n";
             for(my $i=0;$i<=$scaflen{$scaf}-$win;$i=$i+10000){
                   my $start=$i+1;
                   my $end=$i+$win;
                   my $tepercent=sumte($start,$end,$scafte{$scaf});
                   if ($tepercent <= 0.1){
                        print OUT2 "$scaf\t$scaflen{$scaf}\t$tepercent\t$start\t$end\n";
                   } 
             }
               
         }
      }else{
         if ($scaflen{$scaf} < 100000){
             print OUT1 "$scaf\t$scaflen{$scaf}\t0\n";
         }else{
             print OUT1 "$scaf\t$scaflen{$scaf}\t0\n";
             print OUT2 "$scaf\t$scaflen{$scaf}\t0\n";
         }
      }
}
close OUT1;
close OUT2;


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
