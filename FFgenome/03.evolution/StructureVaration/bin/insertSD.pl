#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"psl:s","inf:s","project:s","help");


my $help=<<USAGE;
perl $0 --psl --inf
--inf FPC infromation of BES and clone.
--psl raw blat result of psl format no header
--project: prefix of psl file
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $fpcinf=readfpcinf($opt{inf});
my $title=$1 if ($opt{psl}=~/(.*)\.psl/);
$opt{project}||=$title;
#`perl /share/raid12/chenjinfeng/tools/bin/bestAlign.pl $opt{psl} --cutoff 0.99 > $opt{project}.best.psl`;
`pslReps -singleHit -minCover=0.95 -minAli=0.95 brachyantha.bes2FF.psl brachyantha.bes2FF.filter.psl brachyantha.bes2FF.filter.psr`;

#####stat best hit psl
my $stat=statpsl("$opt{project}.filter.psl");
my ($total,$pe,$pe1,$se);
foreach my $clone (keys %$stat){
    #print "$clone\t$stat->{$clone}->[0]\t$stat->{$clone}->[1]\n";
    $total++;
    if (exists $stat->{$clone}->{0} and exists $stat->{$clone}->{1}){
      $pe++;
      if ($stat->{$clone}->{0}->[0] eq $stat->{$clone}->{1}->[0] and $stat->{$clone}->{0}->[1] ne $stat->{$clone}->{1}->[1] ){
         my $size=$stat->{$clone}->{0}->[3] > $stat->{$clone}->{1}->[2] ? $stat->{$clone}->{0}->[3]-$stat->{$clone}->{1}->[2] : $stat->{$clone}->{1}->[3]-$stat->{$clone}->{0}->[2]; 
         push (@insert,$size);
         #print "$clone\t$stat->{$clone}->{0}->[0]\t$stat->{$clone}->{0}->[2]\t$stat->{$clone}->{0}->[3]\t$stat->{$clone}->{1}->[0]\t$stat->{$clone}->{1}->[2]\t$stat->{$clone}->{1}->[3]\n";
         $pe1++;
      }
    }else{
       $se++;
    }
}

my ($mean,$sd,$number)=ci(\@insert);

####print summary
my $singleBES=$fpcinf->{$opt{project}}->[5]+$fpcinf->{$opt{project}}->[6]-2*$fpcinf->{$opt{project}}->[7];
print "Species\tGenome\tDesinature\tClone with PE\tClone with SE\tClone with PE mapped\tClone with SE mapped\n";
print "$fpcinf->{$opt{project}}->[0]\t$fpcinf->{$opt{project}}->[1]\t$fpcinf->{$opt{project}}->[2]\t$fpcinf->{$opt{project}}->[7]\t$singleBES\t$pe\t$se\n";

print "Final PE:$pe1\tMean Size:$mean\tSD:$sd\n";


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

sub statpsl
{
my ($file)=@_;
my %clone;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $name=$unit[9];
    if ($name=~/\w+\_\_\w+?(a\d+\D+\d+)\.(\w+)/){
       my $clone=$1;
       my $direction=$2;
       #print "$clone\t$direction\n";
       if ($direction eq "r"){
           $clone{$clone}{1}=[$unit[13],$unit[8],$unit[15],$unit[16]];
       }elsif($direction eq "f"){
           $clone{$clone}{0}=[$unit[13],$unit[8],$unit[15],$unit[16]];
       }
    }
}
close IN;
return \%clone;
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

