#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"list:s","type:s","col:s","help");


my $help=<<USAGE;
Summary gene type for gene classified in input list.
perl $0 --list --type
--list: 
--type: 
--col: colume of which the list will be summarized
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{col} ||= 1;
my $reftype=readtable($opt{type});
my $reflist=readtable2($opt{list});

my %hash;
foreach my $g (keys %$reflist){
   my %temp=(
    "Total" => 0,
    "notfound" => 0,
    "noncluster" => 0,
    "nonortholog" => 0,
    "nonsyn.ortholog" => 0,
    "syn.ortholog" => 0,
    "unique" => 0,
    "tandem" => 0
   );
   $hash{$reflist->{$g}}=\%temp unless (exists $hash{$reflist->{$g}});
   $hash{$reflist->{$g}}{"Total"}++;
   if (exists $reftype->{$g}){
      $hash{$reflist->{$g}}{$reftype->{$g}}++;
   }else{
      $hash{$reflist->{$g}}{"notfound"}++;
   }
}

  my %temp=(
    "Total" => 0,
    "notfound" => 0,
    "noncluster" => 0,
    "nonortholog" => 0,
    "nonsyn.ortholog" => 0,
    "syn.ortholog" => 0,
    "unique" => 0,
    "tandem" => 0
   );

my @temp1=sort {$a cmp $b} keys %temp;
my $head=join("\t",@temp1);
print "Type\t$head\n";

foreach my $c (sort keys %hash){
   print "$c\t";
   my @temp=sort {$a cmp $b} keys %{$hash{$c}};
   foreach my $t (@temp){
      print "$hash{$c}->{$t}\t";
   }
   print "\n";
}


########################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if $unit[0]=~/LOC_(.*)\.\d+/;
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}


sub readtable2
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if $unit[0]=~/LOC_(.*)\.\d+/;
    $hash{$unit[0]}=$unit[$opt{col}];
}
close IN;
return \%hash;
}

