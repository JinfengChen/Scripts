#!/usr/bin/perl

if (@ARGV < 1){
   print "Please input a chr.scaffold file: perl mergescaf.pl chr.scaffold > log \n";
   exit();
}

##### merge scaffold sequence. gene gff and te gff file into chromosome by chr.scaffold
my $scaffold="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/seq/OBa.all.fa.RepeatMasker.masked";
my $genegff="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OBa.all.gff";
my $tegff="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OBa.all.fa.RepeatMasker.out.gff";
my $manualte="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OBa.all.manual.TE.gff";
my $rootqry="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/Root.bed";
my $shootqry="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/Shoot.bed";
my $methy="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/FF.methylation.bed";
my $h3k4="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OB.H3K4.bed";
my $h3k4quarter="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OB.H3K4.quarter.bed";
my $scafseq=fastaseq($scaffold);
my %chrseq;
my %scaf2chr;

my $chrscaffold=$ARGV[0];
open IN,"$chrscaffold" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if (exists $chrseq{$unit[0]}){
        $chrseq{$unit[0]}.=$scafseq->{$unit[1]};
    }else{
        $chrseq{$unit[0]}=$scafseq->{$unit[1]};
    }
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
##################write linked seqence to file
open OUT, ">/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/seq/ffseq.masked.chr" or die "$!";
foreach (keys %chrseq){
    print OUT ">$_\n$chrseq{$_}\n"
}
close OUT;
undef %chrseq;
undef $scafseq;
###############################################
########################################################################
#merge($tegff,\%scaf2chr);
merge($genegff,\%scaf2chr);
#merge($manualte,\%scaf2chr);
#mergesoap($rootqry,\%scaf2chr);
#mergesoap($shootqry,\%scaf2chr);
#mergesoap($methy,\%scaf2chr);
#mergesoap($h3k4,\%scaf2chr);
#mergesoap($h3k4quarter,\%scaf2chr);
############soap merge##############################################
sub mergesoap {
my ($gff,$scaf2chr)=@_;
my $dir= $gff.".chr";
`mkdir $dir` if (!-d $dir);
open IN, "$gff" or die "$!";
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
          $temp[2]+=$scafstart;
       }else{
          my $end   =$scafstart+$scaflen-$temp[1];
          my $start =$scafstart+$scaflen-$temp[2];
          $temp[1]=$start;
          $temp[2]=$end;
          if ($temp[5] eq "-"){
                   $temp[5] = "+";
          }elsif($temp[5] eq "+"){
                   $temp[5] = "-"; 
          }     
       }
       $temp[0]=$scafref;
       my $line=join("\t",@temp); 
       open OUT, ">>$dir/$scafref" or die "$!"; 
            print OUT "$line\n"; 
       close OUT;
     }
    }else{
       $unit[0]=~s/Scaffold000/Super/;
       my $line=join("\t",@unit);
       open OUT, ">>$dir/chrUN" or die "$!";
            print OUT "$line\n";
       close OUT;
    }
}
close IN;
}

############te merge#####################################################
sub merge {
my ($gff,$scaf2chr)=@_;
my $dir= $gff.".chr";
`mkdir $dir` if (!-d $dir);
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    if ($_ eq ""){
       print OUT "$_\n";
       next;
    }
    if ($_ =~/^##/){
       print OUT "$_\n";
       next;
    }
    my @unit=split("\t",$_);
    my $scaf=$unit[0];
    if (exists $scaf2chr->{$unit[0]}){
     my $array=$scaf2chr->{$unit[0]};
     foreach (@$array){
       my @temp=@unit;
       my $scafstrand=$_->[3];
       my $scafstart =$_->[1];
       my $scaflen   =$_->[2];
       my $scafref   =$_->[0]; 
       #print "$scafref\t$scaf\t$scaflen\t$scafstart\n";
       if ($scafstrand eq "+"){
          $temp[3]+=$scafstart;
          $temp[4]+=$scafstart;
       }else{
          my $end   =$scafstart+$scaflen-$temp[3];
          my $start =$scafstart+$scaflen-$temp[4];
          $temp[3]=$start;
          $temp[4]=$end;
          if ($temp[6] eq "-"){
                   $temp[6] = "+";
          }elsif($temp[6] eq "+"){
                   $temp[6] = "-"; 
          }     
       }
       #print "$scafref\t$scaf\t$scaflen\t$scafstart\t$temp[3]\n";
       $temp[0]=$scafref;
       my $line=join("\t",@temp); 
       open OUT, ">>$dir/$scafref" or die "$!"; 
            print OUT "$line\n"; 
       close OUT;
     }
    }else{
       $unit[0]=~s/Scaffold000/Super/;
       my $line=join("\t",@unit);
       open OUT, ">>$dir/chrUN" or die "$!";
            print OUT "$line\n";
       close OUT;
    }
}
close IN;
}
###############################################################3333
sub fastaseq {
my ($file)=@_;
my %seq;
$/=">";
print "$file\n";
open IN, "$file" or die "$!";
      while (<IN>){
          chomp $_;
          next if (length $_ < 2);
          my @unit=split("\n",$_);
          my $head=shift @unit;
          my $seq=join("",@unit);
          $seq=~s/\s//g;
          $seq=~s/\r//g;
          $seq=~s/\>//g;
          $seq{$head}=$seq;
      }
close IN;
$/="\n";
return \%seq;
}
