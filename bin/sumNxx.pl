## this scripts read in a fasta file of scaffold sequence
## it can split the sequence by NN to generate a contig file
## then calculate the N50 and N90 of both scaffold and contig to output

die "Usage: perl sumNxx.pl infile > summary\n" if (@ARGV < 1);
my $file=$ARGV[0];
$/="\>";
my $totallen;
my $scaftotal;
my @len;
my @scaflen;
open IN, "$file" or die "can not open my infile";
      while (<IN>){ 
            next if (length $_ < 2);
            my @unit=split("\n",$_);
            my $head=shift @unit;
            #print "$head\n";
            my $seq=join("",@unit);
            my $scafl=length $seq;
            #next if ($scafl < 300000);
            $scaftotal+=$scafl;
            push (@scaflen,$scafl);
            my @contig=split(/N+/,$seq);
            foreach (@contig){
                  $counter++;
                  my $l=length $_;
                  $totallen+=$l;
                  push (@len,$l);
                  #print ">contig$counter\n$_\n";
            }
      }
close IN;
$/="\n";
print "Scaffold Summary\n";
show($scaftotal,\@scaflen); ### sum scaffold infor
print "Contig Summary\n";
show($totallen,\@len); ### sum contig infor

sub show {
my ($totallen,$len)=@_;
my $len50=$totallen*0.5;
my $len60=$totallen*0.6;
my $len70=$totallen*0.7;
my $len80=$totallen*0.8;
my $len90=$totallen*0.9;
my $len95=$totallen*0.95;
print "Total Size: $totallen\n";
print "NXX: Size\tSumLen\tNumber\n";
my @unit=tellme($len50,$len);
print "N50: $unit[0]\t$unit[1]\t$unit[2]\n";
my @unit=tellme($len60,$len);
print "N60: $unit[0]\t$unit[1]\t$unit[2]\n";
my @unit=tellme($len70,$len);
print "N70: $unit[0]\t$unit[1]\t$unit[2]\n";
my @unit=tellme($len80,$len);
print "N80: $unit[0]\t$unit[1]\t$unit[2]\n";
my @unit=tellme($len90,$len);
print "N90: $unit[0]\t$unit[1]\t$unit[2]\n";
my @unit=tellme($len95,$len);
print "N95: $unit[0]\t$unit[1]\t$unit[2]\n";

my @array=tellme2($len);
print "Longest: $array[0]\t> 100bp: $array[1]\t> 2kb: $array[2]\n";

}

sub tellme2
{
my ($len)=@_;
my @temp=sort {$b <=> $a} @$len;
my $bp;
my $kb;
foreach (@temp){
   if ($_ >= 100){
     $bp++;
   }
   if ($_ >= 2000){
     $kb++;
   }
}
my $longest=shift @temp;
return ($longest,$bp,$kb);
}


sub tellme {
my ($len50,$len)=@_;
my $sum1;
my $n50;
my $counter;
foreach (sort {$b <=> $a} @$len){
      $sum1+=$_;
      $counter++;
      if ($sum1 > $len50){
          #print "N50: $_\n";
          $n50=$_;
          last;
      }
}

return ($n50,$sum1,$counter);

}

=pod
foreach (sort {$b <=> $a} @len){
      $sum2+=$_;
      if ($sum2 > $len90){
          print "N90: $_\n";
          last; 
      }
}
=cut



