#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","help");


my $help=<<USAGE;

perl $0 -fasta ff.fasta > summary
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $total;
my $hypo;
my $similar;
my $other;
my $refanno=getfastaanno("$opt{fasta}");
foreach (keys %$refanno){
     #print "$_\t$refanno->{$_}\n";
     $total++;
     if ($refanno->{$_}=~/hypothetical/i){
     #if ($refanno->{$_}=~/hypothetical/i or $refanno->{$_}=~/predicted/i){
         $hypo++;
         print "Hypo\t$_\t$refanno->{$_}\n";
     }elsif($refanno->{$_}=~/similar/i){
         $similar++;
         print "Similar\t$_\t$refanno->{$_}\n";
     }else{
         $other++;
         print "Other\t$_\t$refanno->{$_}\n";
     }
}

print "Total: $total\n";
print "Hypothetical: $hypo\n";
print "Similar: $similar\n";
print "Other: $other\n";
#####################################
sub getfastaanno
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
    my $head= shift @temp1;
    my $anno= join(" ",@temp1);
    #print "$head\n";
    $hash{$head}=$anno;
}
$/="\n";
return \%hash;
}
 
