#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","link:s","project:s","help");


my $help=<<USAGE;
Generate a chr.scaffold table that used in drawFeature from link and super-scaffold.fa file.
perl $0 --fasta --link
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}
my %hash;
my %chrlen;
my %temp;
my $reflen = getfastalen($opt{fasta});
open IN, "$opt{link}" or die "$!";
<IN>;
while (<IN>){
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $chrname;
    if ($unit[0]=~/chr(\d+)/){
       $unit[0]=$1;
    }else{
       next;
    }
    next if (exists $temp{$unit[1]});
    $temp{$unit[1]}=1;
    my $scafname=sprintf("Scaffold%06d",$unit[1]);
    if (exists $chrlen{$unit[0]}){
       $chrlen{$unit[0]}+=$reflen->{$scafname};
       $chrlen{$unit[0]}+=100;
    }else{
       $chrlen{$unit[0]}+=$reflen->{$scafname};
    }
    my $start=$chrlen{$unit[0]}-$reflen->{$scafname};
    push (@{$hash{$unit[0]}},[$scafname,$start,$chrlen{$unit[0]}]);
}
close IN;

open OUT, ">OBa_chr_scaffold" or die "$!";
print OUT "Chr\tScaffold\tScafLen\tScaffoldStartONChr\tScaffoldEndOnChr\tScaffoldStrand\n";
foreach my $chr (sort {$a <=> $b} keys %hash){
    $chrname=sprintf("chr%02d",$chr);
    foreach my $line (sort {$a->[1] <=> $b->[1]} @{$hash{$chr}}){
        #my $scafname=sprintf("Scaffold%06d",$scaf);
        print OUT "$chrname\t$line->[0]\t$reflen->{$line->[0]}\t$line->[1]\t$line->[2]\t\+\n";
    }
}
close OUT;

####################
sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    my $len=length $seq;
    $hash{$head}=$len;
}
$/="\n";
return \%hash;
}

