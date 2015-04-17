#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
our %opt;
GetOptions (\%opt,"ortholog:s","tandem:s","transpose1:s","transpose2:s","help");


my $help=<<USAGE;
perl $0 --ortholog --tandem --transpose1 --transpose2

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my (%total1,%total2,%ortholog1,%ortholog2,%tandem1,%tandem2,%transpose1,%transpose2);
my ($orth1,$orth2)=orth($opt{ortholog});
my ($tan1,$tan2)=tandem($opt{tandem});
my $trans1=transpose($opt{transpose1});
my $trans2=transpose($opt{transpose2});
open OUT, ">manaul.OBa.RLKLRR.ID" or die "$!";
foreach my $g (sort keys %$orth2){
   #print "$g\n";
   $total2{$orth2->{$g}}++;
   $ortholog2{$orth2->{$g}}++;
   print OUT "$g\t$orth2->{$g}\n";
}
foreach my $g (sort keys %$tan2){
   #print "$g\n";
   $total2{$tan2->{$g}}++;
   $tandem2{$tan2->{$g}}++;
   print OUT "$g\t$tan2->{$g}\n";
}
foreach my $g (sort keys %$trans2){
   $total2{$trans2->{$g}}++;
   $transpose2{$trans2->{$g}}++;
   print OUT "$g\t$trans2->{$g}\n";
}
close OUT;
 
open OUT, ">manaul.rice.RLKLRR.ID" or die "$!";
foreach my $g (sort keys %$orth1){
   $total1{$orth1->{$g}}++;
   $ortholog1{$orth1->{$g}}++;
   print OUT "$g\t$orth1->{$g}\n";
}
foreach my $g (sort keys %$tan1){
   $total1{$tan1->{$g}}++;
   $tandem1{$tan1->{$g}}++;
   print OUT "$g\t$tan1->{$g}\n";
}
foreach my $g (sort keys %$trans1){
   $total1{$trans1->{$g}}++;
   $transpose1{$trans1->{$g}}++;
   print OUT "$g\t$trans1->{$g}\n";
}
close OUT;
my @total;
print "Class\tRice\tOrtholog\tTandem\tTranspose\tGramene\tOrtholog\tTandem\tTranspose\n";
foreach my $type (sort keys %ortholog1){
   $tandem1{$type}= $tandem1{$type} ? $tandem1{$type} : "0";
   $tandem2{$type}= $tandem2{$type} ? $tandem2{$type} : "0";
   $transpose1{$type} = $transpose1{$type} ? $transpose1{$type} : "0";
   $transpose2{$type} = $transpose2{$type} ? $transpose2{$type} : "0";
   print "$type\t$total1{$type}\t$ortholog1{$type}\t$tandem1{$type}\t$transpose1{$type}\t$total2{$type}\t$ortholog2{$type}\t$tandem2{$type}\t$transpose2{$type}\n";
}

##################
sub orth
{
my ($file)=@_;
my %hash1;
my %hash2;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print "$unit[6]\n" if exists $hash{$unit[6]};
    $hash2{$unit[6]}=$unit[1];
    $hash1{$unit[0]}=$unit[1];
}
close IN;
return (\%hash1,\%hash2);
}

sub transpose
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub tandem
{
my ($file)=@_;
my %hash1;
my %hash2;
my $tandem1;
my $tandem2;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=~s/\"//g;
    my @ref=split(",",$unit[0]);
    $unit[6]=~s/\"//g;
    my @qry=split(",",$unit[6]);
    #print "$unit[0]\n$unit[6]\n";
    foreach my $g (@qry){
       $hash2{$g}=$unit[1];
    }
    foreach my $g (@ref){
       $hash1{$g}=$unit[1];
    }
    $tandem1+=@ref;
    $tandem2+=@qry;
}
close IN;
print "Ref tandem:$tandem1\nQry tandem:$tandem2\n";
return (\%hash1,\%hash2);
}
       

