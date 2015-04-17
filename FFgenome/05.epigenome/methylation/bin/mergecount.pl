#!/usr/bin/perl

if (@ARGV < 1){
   print "Please input a chr.scaffold file: perl $0 chr.scaffold > log \n";
   exit();
}

##### merge scaffold sequence. gene gff and te gff file into chromosome by chr.scaffold
my $methylation="../input/FF.cout";
my %scaf2chr;

my $chrscaffold=$ARGV[0];
open IN,"$chrscaffold" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if (exists $scaf2chr{$unit[1]}){
        #print "$unit[1]\t$scaf2chr{$unit[1]}\t$unit[0]\n";
        my $temp=$scaf2chr{$unit[1]};
        push (@$temp,[$unit[0],$unit[3],$unit[2],$unit[5]]);
        $scaf2chr{$unit[1]}=$temp;
    }else{
        my @temp;
        push (@temp,[$unit[0],$unit[3],$unit[2],$unit[5]]);
        $scaf2chr{$unit[1]}=\@temp;
    }
}
close IN;
########################################################################
mergecout($methylation,\%scaf2chr);
############soap merge##############################################
sub mergecout {
my ($cout,$scaf2chr)=@_;
my $dir= $cout.".chr";
`mkdir $dir` if (!-d $dir);
open IN, "$cout" or die "$!";
while (<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if (exists $scaf2chr->{$unit[0]}){
     my $array=$scaf2chr->{$unit[0]};
     foreach (@$array){
       my @temp=@unit;
       my $scafstrand=$_->[3];
       my $scafstart =$_->[1];
       my $scaflen   =$_->[2];
       my $scafref   =$_->[0]; 
       if ($scafstrand eq "+"){
          $temp[1]+=$scafstart;
       }else{
          my $pos   =$scafstart+$scaflen-$temp[1];
          $temp[1]=$pos;
          if ($temp[5] eq "-"){
                   $temp[5] = "+";
          }elsif($temp[5] eq "+"){
                   $temp[5] = "-"; 
          }     
       }
       $temp[0]=$scafref;
       my $line=join("\t",@temp); 
       open OUT, ">>$dir/$scafref.cout" or die "$!"; 
            print OUT "$line\n"; 
       close OUT;
     }
    }else{
       $unit[0]="chrUN";
       my $line=join("\t",@unit);
       open OUT, ">>$dir/chrUN.cout" or die "$!";
            print OUT "$line\n";
       close OUT;
    }
}
close IN;
}

