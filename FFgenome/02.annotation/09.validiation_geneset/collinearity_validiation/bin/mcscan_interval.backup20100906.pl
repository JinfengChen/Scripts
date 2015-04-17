#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"refgff:s","refpep:s","qrygff:s","qrypep:s","align:s","help");


my $help=<<USAGE;

perl $0 -refgff Os.gff -refpep Os.pep -qrygff Ob.gff -qrypep Ob.pep -align ob_os.align 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


$blastall="/home/biosoftware/blast-2.2.14/bin/blastall";
$formatdb="/home/biosoftware/blast-2.2.14/bin/formatdb";
$bestalign="/home/jfchen/FFproject/tools/bin/bestAlign.pl";



my $refgff=gff($opt{refgff});
our $refpep=getfastaseq($opt{refpep});
my $qrygff=gff($opt{qrygff});
our $qrypep=getfastaseq($opt{qrypep});
my $align =parsealign($opt{align});

my $additionalpair;

open INTER, ">interval.txt" or die "$!";
foreach(sort {$a <=> $b} keys %$align){
	 my $block=$_; 
	 print "Block: $block\n";
	 my $reftemp=$align->{$_};
         my $chr=shift @$reftemp;
         my $qrychrgff=$qrygff->{$chr->[0]};
         my $refchrgff=$refgff->{$chr->[1]};
         print "$chr->[0]\t$chr->[1]\n";
         for(my $i=1;$i<@$reftemp;$i++){
             $counter++;
             my $last=$reftemp->[$i-1];
             my $present=$reftemp->[$i];
             print "Pair: $counter\n";
             print "Last: $last->[0]\t$last->[1]\n";
             print "Present: $present->[0]\t$present->[1]\n";
             ####deal with query####
             my $qrystart=$qrychrgff->{$last->[0]};
	     my $qryend  =$qrychrgff->{$present->[0]};
	     my $qrygene =countgene($qrychrgff,$qrystart,$qryend);
             my $tempn1=@$qrygene;
             print "$tempn1\n@$qrygene\n";
             print INTER "$tempn1\t";
             ####deal with refer#### 
             my $refstart=$refchrgff->{$last->[1]};
             my $refend  =$refchrgff->{$present->[1]};
             my $refgene =countgene($refchrgff,$refstart,$refend);
             my $tempn2=@$refgene;
             print "$tempn2\n@$refgene\n";
             print INTER "$tempn2\n";
             ####ortholog pair
             #my $orthpair=ortholog ($qrygene,$refgene);
             #foreach (keys %$orthpair){
             #      print "$_\t$orthpair->{$_}\n";
             #}
             my $temp=dealtwolist ($qrygene,$refgene);
             $additionalpair+=$temp;
         }
}
close INTER;

print "Get additional pair of ortholog: $additionalpair\n";

###############sub functions####################

sub dealtwolist
{
my ($qrygene,$refgene)=@_;
my $additionalpair=-2;
& getlistseq ($qrygene,$qrypep,"qry.pep");
& getlistseq ($refgene,$refpep,"ref.pep");
system ("$formatdb -i qry.pep -p T");
system ("$formatdb -i ref.pep -p T");
system ("$blastall -p blastp -i qry.pep -d qry.pep -e 1e-5 -o qry_qry.blastm8 -m 8 > qry.log 2> qry.log2");
system ("$blastall -p blastp -i ref.pep -d ref.pep -e 1e-5 -o ref_ref.blastm8 -m 8 > ref.log 2> ref.log2");
my $qrymost=mosthitm8("qry_qry.blastm8");
my $refmost=mosthitm8("ref_ref.blastm8");
if ($qrymost >= 2 or $refmost >= 2){
   print "Have tandem duplicate\n";
   `rm qry* ref*`;
}else{
   my $orthpair=ortholog();
   foreach (keys %$orthpair){
      $additionalpair++;
      print "$_\t$orthpair->{$_}\n";
   }
}
return $additionalpair;
}


sub ortholog
{
my %best_pair1;
my %best_pair2;
my %best_pair;
system ("$blastall -p blastp -i qry.pep -d ref.pep -e 1e-5 -o qry_ref.blastm8 -m 8 > qry.log 2> qry.log2");
system ("$blastall -p blastp -i ref.pep -d qry.pep -e 1e-5 -o ref_qry.blastm8 -m 8 > ref.log 2> ref.log2");
system ("perl $bestalign -f m8 qry_ref.blastm8 > qry_ref.blastm8.best");
system ("perl $bestalign -f m8 ref_qry.blastm8 > ref_qry.blastm8.best");
& read_m8_best_pair ("qry_ref.blastm8.best",\%best_pair1);
& read_m8_best_pair ("ref_qry.blastm8.best",\%best_pair2);
`rm qry* ref*`;
foreach my $pep1 (sort keys %best_pair1) {
        my $pep2 = $best_pair1{$pep1};
        if (exists $best_pair2{$pep2} && $best_pair2{$pep2} eq $pep1) {
           $best_pair{$pep1} = $pep2;
        }
}
return \%best_pair;
}


sub mosthitm8
{
my ($m8)=@_;
my $most=0;
my %hash;
open IN, "$m8" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   if (exists $hash{$unit[0]}){
       $hash{$unit[0]}++;
   }else{
       $hash{$unit[0]}=0; ### at least have a self hit
   }
}
close IN;
foreach(keys %hash){
   if ($hash{$_} > $most){
      $most=$hash{$_};
   }
}
return $most;
}

sub getlistseq
{
my ($refarray,$refseq,$outfile)=@_;
open TEMP, ">$outfile" or die "$!";
foreach(@$refarray){
    if (exists $refseq->{$_}){
        print TEMP ">$_\n$refseq->{$_}\n";
    }
}
close TEMP;
}

sub read_m8_best_pair {
        my $file = shift;
        my $hash_p = shift;
        open IN,$file || die "fail $file";
        while (<IN>) {
                my @t = split /\t/;
                $hash_p->{$t[0]} = $t[1];
        }
        close IN;
}



sub getfastaseq
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
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


sub countgene
{
my ($chrgff,$start,$end)=@_;
if ($start > $end){
   my $temp=$start;
   $start= $end;
   $end  = $temp;
}
print "$start\t$end\n";
my $number=0;
my @gene=sort {$chrgff->{$a} <=> $chrgff->{$b}} keys %$chrgff;
my @count;
foreach(@gene){
    if ($chrgff->{$_} >= $start and $chrgff->{$_} <= $end) {
        $number++;
        push (@count,$_);
    }elsif($chrgff->{$_} > $end){
	return \@count;
    }
}
return \@count;
}





sub parsealign
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
	   push (@$reftemp,[$unit[1],$unit[2]]);
	   $hash{$block}=$reftemp;		
	}else{
           my @temp;
           my @unit=split("\t",$_);
           push (@temp,[$qry,$ref]);
	   push (@temp,[$unit[1],$unit[2]]);
           $hash{$block}=\@temp;
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
