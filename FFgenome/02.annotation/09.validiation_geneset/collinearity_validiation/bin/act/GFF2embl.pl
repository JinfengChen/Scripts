#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","embl:s","fasta:s","rev","help");


my $help=<<USAGE;
convert format between gff and embl.
-rev: if this option is used, convert embl to gff
-fasta: if this option is used, sequence will be add to embl when converting gff to embl
perl GFF2embl.pl -gff OB.gff -embl OB.embl -rev
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

& gff2embl($opt{gff},$opt{embl});


#####################
sub gff2embl
{
my ($gff,$embl)=@_;
my $refgff=parseGFF($gff);
open OUT, ">$embl" or die "$!";
print OUT "FH   Key             Location\/Qualifiers\n";
foreach(sort keys %$refgff){
    my $id=$_;
    my @line=split("\n",$refgff->{$id});
    my $mRNA=shift @line;
    my @mRNA=split("\t",$mRNA);
    my $strand=$mRNA[6];
    my @cds;
    foreach(@line){
       my @unit=split("\t",$_);
       push (@cds,[$unit[3],$unit[4]]);
    }
    my @sortcds=sort {$a->[0] <=> $b->[0]} @cds;
    my @cdspos;
    foreach (@sortcds){
       my $temp=$_->[0]."..".$_->[1];
       push (@cdspos,$temp);
    }
    my $pos=join(",",@cdspos);
    #print "$strand\n";
    if ($strand eq "+"){
       print OUT "FT   mRNA            join($pos)\n";
       print OUT "FT                   \/gene\=\"$id\"\n";
       print OUT "FT   CDS             join($pos)\n";
       print OUT "FT                   \/gene\=\"$id\"\n";
    }else{
       print OUT "FT   mRNA            complement(join($pos))\n";
       print OUT "FT                   \/gene\=\"$id\"\n";
       print OUT "FT   CDS             complement(join($pos))\n";
       print OUT "FT                   \/gene\=\"$id\"\n";
    }
}
if (exists $opt{fasta}){
        my $refseq=getfastaseq($opt{fasta});
        my @head=keys %$refseq;
        my $seq=$refseq->{$head[0]};
        $len = length ($seq);
        $mod = $len % 60;
        $times = ($len - $mod) / 60;
        print OUT "SQ   Sequence $len BP;\n";
        for ($i = 0; $i < $times; $i++) {
                $a = substr ($seq, $i * 60, 60);
                print OUT "    ";
                for ($j = 0; $j < 6; $j++) {
                        print OUT " ", substr ($a, $j * 10, 10);
                }
                printf OUT ("%10d\n", $i * 60 + 60);
        }
        $a = substr ($seq, $len - $mod, $mod);
        print OUT "    ";
        for ($j = 0; $j < 6; $j++) {
                print OUT " ", substr ($a, $j * 10, 10);
                $l = $l." ".substr ($a, $j * 10, 10);
        }
        $l = length ($l) + 4;
        $l = 70 - $l;
        for ($j = 0; $j < $l; $j++) {
                print OUT " ";
        }
        printf OUT ("%10d\n", $i * 60 + $mod);
}
print OUT "\/\/\n";
close OUT;
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
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
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


