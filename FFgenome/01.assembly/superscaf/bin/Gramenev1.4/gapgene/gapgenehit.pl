#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"pep:s","gff:s","blast:s","list:s","help");


my $help=<<USAGE;
Check blasttable and summary the hit information for two part of gap gene.
perl $0 --gff --blast --list
--pep: pep sequence of two part of gap gene
--gff: gff file for print out the position of gene
--blast: blasttable
--list: check list of 100 gap gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $gff=parseGFF($opt{gff});
my $list=getlist($opt{list});
my $length=getfastalen($opt{pep});
my ($up,$down)=getblast($opt{blast});

print "rank\tgene\tgap length\tposition\tcoverage\tfive\tup hit length\tup hit\tup subject length\tup identity\tthree\tdown hit length\tdown hit\tdown subject length \tdown identity\n";
foreach my $rank (keys %$list){
    foreach my $gene (sort keys %{$list->{$rank}}){
        my $gaplen=$list->{$rank}->{$gene};
        my $five=$length->{$gene}{5} ? $length->{$gene}{5} : 0;
        my $three=$length->{$gene}{3} ? $length->{$gene}{3} : 0;
        my $position=$gff->{$gene};
        my ($uphitlen,$uphit,$upsubjectlen,$upid);
        if (exists $up->{$gene}){
           $uphitlen=$up->{$gene}->[0];
           $uphit=$up->{$gene}->[1];
           $upsubjectlen=$up->{$gene}->[2];
           $upid=$up->{$gene}->[3];
        }else{
           $uphitlen=0;
           $uphit="NA";
           $upsubjectlen="NA";
           $upid="NA";
        }
        my ($downhitlen,$downhit,$downsubjectlen,$downid);
        if (exists $down->{$gene}){
           $downhitlen=$down->{$gene}->[0];
           $downhit=$down->{$gene}->[1];
           $downsubjectlen=$down->{$gene}->[2];
           $downid=$down->{$gene}->[3];
        }else{
           $downhitlen=0;
           $downhit="NA";
           $downsubjectlen="NA";
           $downid="NA";
        }
        my $subjectlen=0;
        if ($upsubjectlen ne "NA"){
           $subjectlen=$upsubjectlen;
        }
        if ($downsubjectlen ne "NA"){
           $subjectlen=$downsubjectlen;
        }
        my $coverage=$subjectlen > 0 ? ($uphitlen+$downhitlen)/$subjectlen : 0;
        $coverage= $coverage > 1 ? 1 : $coverage;
        print "$rank\t$gene\t$gaplen\t$position\t$coverage\t$five\t$uphitlen\t$uphit\t$upsubjectlen\t$upid\t$three\t$downhitlen\t$downhit\t$downsubjectlen\t$downid\n";
    }
}

#################################################

sub getblast{
my ($file)=@_;
my (%up,%down);
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[0]=~/(.*)\_3/){
       my $hitlen=$unit[3]-$unit[2]+1;
       my $hit=$unit[4];
       my $identity=$unit[8];
       my $subjectlen=$unit[5];
       $down{$1}=[$hitlen,$hit,$subjectlen,$identity];
    }elsif($unit[0]=~/(.*)\_5/){
       my $hitlen=$unit[3]-$unit[2]+1;
       my $hit=$unit[4];
       my $identity=$unit[8];
       my $subjectlen=$unit[5];
       $up{$1}=[$hitlen,$hit,$subjectlen,$identity];
    }
}
close IN;
return (\%up,\%down);
}


sub getlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $hash{$unit[0]}{$unit[1]}=$unit[2]; ## 1->gene->gaplen or 2->gene->gaplen
}
close IN;
return \%hash;
}

sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my ($gene,$part);
    if ($head=~/(.*)\_(\d+)/){
       $gene=$1;
       $part=$2;
    }
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$gene}{$part}=length $seq; # gene->5->len or gene->3->len
}
$/="\n";
return \%hash;
}

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        my $id;
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
            $hash{$id}=$unit[3]; # gene->start
        }
    }

}
close IN;
return \%hash;
}

