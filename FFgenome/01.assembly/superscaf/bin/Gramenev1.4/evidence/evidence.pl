#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gene:s","protein:s","transcript:s","iprscan:s","help");


my $help=<<USAGE;
perl $0 --gene --protein --transcript --iprscan 
--gene: cds fasta
--protein: blasttable file produced by bgi blastparse.pl, top hit 1, top match 1. with header.
--transcript: 
--iprscan: iprscan raw table
Get evidence code for each gene in gramenev1.4. 
code: 1:homolog protein or have Pfam annotation
code: 2:only have RNA seq transcript evidence
code: 3:have no evidence
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $gene=getfastaseq($opt{gene});
my $protein=getblast($opt{protein},0,0);
my $transcript=getblast($opt{transcript},0,0.95);
my $pfam=getiprscan($opt{iprscan});

my %hash;
open OUT, ">evidence.code" or die "$!";
foreach my $g (sort keys %$gene){
    if (exists $protein->{$g} or exists $pfam->{$g}){
        $hash{"type1"}++;
        if (exists $protein->{$g}){
           print OUT "$g\t$gene->{$g}\t1\t$protein->{$g}";
           if ($pfam->{$g}){
              print OUT "\t$pfam->{$g}\n";
           }else{
              print OUT "\n";
           }
        }else{
           print OUT "$g\t$gene->{$g}\t1\tNA\t$pfam->{$g}\n";

        }
    }elsif(exists $transcript->{$g}){
        $hash{"type2"}++;
        print OUT "$g\t$gene->{$g}\t2\t$transcript->{$g}\n"; 
    }else{
        $hash{"type3"}++;
        print OUT "$g\t$gene->{$g}\t3\tNo Evidence\n";
    }
}
close OUT;
foreach (sort keys %hash){
    print "$_\t$hash{$_}\n";
}
#########################################################
sub getiprscan
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
     while(<IN>){
            chomp $_;
            my @unit=split("\t",$_);
            next unless ($unit[3] eq "HMMPfam");
            unless (exists $hash{$unit[0]}){
                $hash{$unit[0]}="$unit[0]\t$unit[12]";
            }
     }
close IN;
return \%hash;
}


sub getblast
{
my ($file,$cov,$id)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[13] > 1e-10);
    if ($unit[0]=~/\w+/){
       my $hitlen=$unit[3]-$unit[2]+1;
       my $hit=$unit[4];
       my $identity=$unit[8];
       my $subjectlen=$unit[5];
       next if ($hitlen/$unit[1] < $cov or $identity < $id);
       $hash{$unit[0]}="$hit:$hitlen:$identity" unless exists $hash{$unit[0]};
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
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}

