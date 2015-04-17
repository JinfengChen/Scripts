#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"refgff:s","qrygff:s","align:s","chrlen:s","help");


my $help=<<USAGE;
Calculate the sequence and gene coverage for each chromosome in Mcscan block.
perl $0 -refgff Os.gff -qrygff Ob.gff -align ob_os.align -chrlen ob_os.chr 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff=gff($opt{refgff});
my $qrygff=gff($opt{qrygff});
my $align=align($opt{align});
my $chrlen=chrlen($opt{chrlen});

my %lencoverage;
my %genecoverage;
my %addlen;
my %addgene;
my %chrgene;
my %cogene; ### colinear gene
### count gene number for each chromosome
foreach(keys %$refgff){
        my $chr=$_;
	my $tempgff=$refgff->{$chr};
	my @gene=keys %$tempgff;
	my $tempnum=@gene;
	print "$chr\t$tempnum\n";
	$chrgene{$chr}=$tempnum;
}

foreach(keys %$qrygff){
        my $chr=$_;
	my $tempgff=$qrygff->{$chr};
	my @gene=keys %$tempgff;
	my $tempnum=@gene;
	print "$chr\t$tempnum\n";
	$chrgene{$chr}=$tempnum;
}
##########################################


=pod
foreach(sort keys %$chrlen){
   print "$_\t$chrlen->{$_}\n";

}
=cut

foreach(sort {$a <=> $b} keys %$align){
	 my $block=$_; 
	 print "Block: $block\n";
	 my $reftemp=$align->{$_};
         foreach(sort {$a cmp $b} keys %$reftemp){
	     my $chr=$_;
             my $colinegene=@{$reftemp->{$chr}};
             print "$chr\t$colinegene\n";
	     my $gene=$reftemp->{$chr};
	     my $startgene=$gene->[0];
	     my $endgene=$gene->[@$gene-1];
	     #print "$chr\t$startgene\t$endgene\n";
	     if (exists $refgff->{$chr}) {
                 my $chrgff=$refgff->{$chr};
	         my $start=$chrgff->{$startgene};
		 my $end  =$chrgff->{$endgene};
                 print "$chr\t$startgene\t$start\t$endgene\t$end\n";
		 my $tempgene=countgene($chr,$chrgff,$start,$end);
                 if (exists $addgene{$chr}) {
		    $addgene{$chr}+=$tempgene;
                 }else{
		    $addgene{$chr}=$tempgene;
		 }
                 if (exists $cogene{$chr}) {
                    $cogene{$chr}+=$colinegene;
                 }else{
                    $cogene{$chr}=$colinegene;
                 }
		 my $templen=abs($end-$start);
                 if (exists $addlen{$chr}) {
		    $addlen{$chr}+=$templen;
                 }else{
	   	    $addlen{$chr}=$templen;
	         }
	     }elsif(exists $qrygff->{$chr}){
		 my $chrgff=$qrygff->{$chr};
	         my $start=$chrgff->{$startgene};
		 my $end  =$chrgff->{$endgene};
		 print "$chr\t$startgene\t$start\t$endgene\t$end\n";
		 my $tempgene=countgene($chr,$chrgff,$start,$end);
                 if (exists $addgene{$chr}) {
	            $addgene{$chr}+=$tempgene;
                 }else{
	            $addgene{$chr}=$tempgene;
	         }
                 if (exists $cogene{$chr}) {
                    $cogene{$chr}+=$colinegene;
                 }else{
                    $cogene{$chr}=$colinegene;
                 }
	         my $templen=abs($end-$start);
                 if (exists $addlen{$chr}) {
		    $addlen{$chr}+=$templen;
                 }else{
		    $addlen{$chr}=$templen;
	         }
	     }
         #print Dumper($gene);
	 }
}


print "Chr\tLength\tCover Length\tCoverage\tGene Number\tCover Number\tCoverage\n";
foreach(sort keys %addlen){
    my $chr=$_;
    my $lencov=$addlen{$chr}/$chrlen->{$chr};
    print "$chr\t$chrlen->{$chr}\t$addlen{$chr}\t$lencov\t";
    my $genecov=$addgene{$chr}/$chrgene{$chr};
    print "$chrgene{$chr}\t$addgene{$chr}\t$genecov\t";
    my $cocov=$cogene{$chr}/$chrgene{$chr};
    print "$cogene{$chr}\t$cocov\n";
}


=pod
foreach(keys %$refgff){
	 my $chr=$_;
     my $chrgff=$refgff->{$_};
	 foreach(keys %$chrgff){
	      print "$chr\t$_\t$chrgff->{$_}\n";
	 }
}
=cut

###############sub functions####################


sub countgene
{
my ($chr,$chrgff,$start,$end)=@_;
my $number=0;
my @gene=sort {$chrgff->{$a} <=> $chrgff->{$b}} keys %$chrgff;
my $chrgene=@gene;
foreach(@gene){
    if ($chrgff->{$_} >= $start and $chrgff->{$_} <= $end) {
        $number++;
    }elsif($chrgff->{$_} > $end){
	return $number;
    }
}
return $number;
}


sub chrlen
{
my ($chr)=@_;
my %hash;
open IN, "$chr" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_!~/^chr/ or $_ eq "");
    my @unit=split("\t",$_);
	$unit[0]=~s/chr/Ob/;
	$unit[2]=~s/chr/Os/;
	$hash{$unit[0]}=$unit[1];
	$hash{$unit[2]}=$unit[3];
}
close IN;
return \%hash;
}



sub align
{
my ($align)=@_;
my %hash;
my $block;
my $qry;
my $ref;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    if ($_=~/\#\# Alignment (\d+)\: score=.* e_value=.* N=\d+ (\w+)\&(\w+) (\w+)/) {
	$block=$1;
	$qry=$2;
	$ref=$3; 
	#print "$block\t$qry\t$ref\n";
    }elsif($_=~/\s*$block\-\s*\d+\:/){
	my @unit=split("\t",$_);
	#print "1:$unit[1]\t2:$unit[2]\n";
	if (exists $hash{$block}) {
           my $reftemp=$hash{$block};
	   push (@{$reftemp->{$qry}},$unit[1]);
           push (@{$reftemp->{$ref}},$unit[2]);
	   $hash{$block}=$reftemp;			
	}else{
           my %temp;
           my @unit=split("\t",$_);
	   push (@{$temp{$qry}},$unit[1]);
           push (@{$temp{$ref}},$unit[2]);
           $hash{$block}=\%temp;
        }
     }	
}
close IN;
return \%hash;
}


sub gff
{
my ($gff)=@_;
my %hash;
open IN, "$gff" or die "$!";
while(<IN>){
        chomp $_;
	next if ($_ eq "");
	my @unit=split("\t",$_);
	if (exists $hash{$unit[0]}) {
        my $reftemp=$hash{$unit[0]};
	    $reftemp->{$unit[1]}=$unit[2];
            $hash{$unit[0]}=$reftemp;
	}else{
	    my %temp;
	    $temp{$unit[1]}=$unit[2];
	    $hash{$unit[0]}=\%temp;
	}
}
close IN;
return \%hash;
}
