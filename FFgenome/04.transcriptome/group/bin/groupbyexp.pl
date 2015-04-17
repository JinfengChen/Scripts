#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"rpkm:s","project:s","help");


my $help=<<USAGE;
Group gene into quintiles by expression level of RPKM (cutoff=1).
perl $0 -rpkm FF.rpkm.tophat -project FF
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my %rank=(
    0 => "1st",
    1 => "2nd",
    2 => "3d",
    3 => "4th",
    4 => "5th"
);


my $refrpkm=rpkmexpr($opt{rpkm});
my @gene=sort {$refrpkm->{$a} <=> $refrpkm->{$b}} keys %$refrpkm;
my $number=int (@gene/5);
my $counter=0;
my %group;
open OUT, ">$opt{project}.gene.group" or die "$!";
foreach (@gene){
   #print "$_\t$refrpkm->{$_}\n";
   $counter++;
   my $index=int ($counter/$number);
   if ($index > 4){
      $index=4;
   }
   if (exists $group{$index}){
       $group{$index}.="\t$_";
   }else{
       $group{$index}=$_;
   }
}
my @temp=sort {$a <=> $b} keys %group;
foreach (@temp){
   print OUT "$rank{$_}\t$group{$_}\n";
}
close OUT;

#############

sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    unless ($unit[1]=~/\-/ or $unit[1] < 1){
       $hash{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}
 
