#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blasttable:s","identity:s","coverage:s","inf:s","project:s","help");


my $help=<<USAGE;
perl $0 --blasttable --inf
--inf FPC infromation of BES and clone.
--blasttable blasttable
--project: 
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

#my $fpcinf=readfpcinf($opt{inf});
my $title=$1 if ($opt{blasttable}=~/(.*)\.blasttable/);
$opt{project}||=$title;

$opt{coverage} ||=0.3;
$opt{identity} ||=0.7;

####
my ($bes,$hit)=stattable("$opt{blasttable}",$title);
my $n1=keys %{$bes->{0}};
my $n2=keys %{$bes->{1}};
print "Total BES with Hit: $n1\n";
print "Total BES Hit: $hit->[0]\n";
print "Total BES with Hit After Filter: $n2\n";
print "Total BES Hit After Filter: $hit->[1]\n";


########################################
sub readfpcinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $name=$1 if ($unit[0]=~/O\. (\w+)/);
    $hash{$name}=[@unit];
}
close IN;
return \%hash;
}

sub stattable
{
my ($file,$title)=@_;
my %bes;
my @hit;
open OUT, ">$title.filter.blasttable" or die "$!";
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    $hit[0]++;
    my @unit=split("\t",$_);
    $bes{0}{$unit[0]}=1;
    $rate=$unit[11]/$unit[1];
    next if ($unit[8] > $opt{identity} and $rate < $opt{coverage});
    print OUT "$_\n";
    $hit[1]++;
    $bes{1}{$unit[0]}=1
}
close IN;
close OUT;
return (\%bes,\@hit);
}


sub ci
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        #print "$_\t$temp\n";
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}
my $number=$loop;
return (0,0,0) if ($number < 2);
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($mean,$SD,$number);
}

