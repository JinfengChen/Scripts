#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"refgff:s","qrygff:s","align:s","help");


my $help=<<USAGE;

perl $0 -refgff Os.gff -qrygff Ob.gff -align ob_os.align 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff=gff($opt{refgff});
my $qrygff=gff($opt{qrygff});
my $align =parsealign($opt{align});


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
             my $last=$reftemp->[$i-1];
             my $present=$reftemp->[$i];
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
         }
}
close INTER;

###############sub functions####################


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
