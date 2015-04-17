#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"block:s","help");


my $help=<<USAGE;
Summary collinearity from block file
perl $0 --block
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$chr= $1 if ($opt{block}=~/(.*)\.blocks/);

my ($refblast,$refgene)=readblast("$chr.blast");
my $refbed=readbed("$chr.bed",$refgene);
& parseblock($opt{block},$refbed,$refblast);
#########

### read blast file, return to hash. one store all genes with blast hit, the other with best hit in self genome.
### hash{gene}->[gene1][e-value]
sub readblast
{
my ($file)=@_;
my %hash;
my %gene;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $gene{$unit[0]}=1;
    $gene{$unit[1]}=1;
    next unless (substr($unit[0],0,2) eq substr($unit[1],0,2));
    #print "$unit[0]\t$unit[1]\n";
    if (exists $hash{$unit[0]}){
       $hash{$unit[0]}=[$unit[1],$unit[2]] if $unit[2] < $hash{$unit[0]}->[1];
    }else{
       $hash{$unit[0]}=[$unit[1],$unit[2]];
    }
    if (exists $hash{$unit[1]}){
       $hash{$unit[1]}=[$unit[0],$unit[2]] if $unit[2] < $hash{$unit[1]}->[1];
    }else{
       $hash{$unit[1]}=[$unit[0],$unit[2]];
    }
}
close IN;
return (\%hash,\%gene);
}

### read bed file, store only these gene with blast hit in hash. hash{gene}->[chr][rank]
sub readbed
{
my ($file,$gene)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    next if (!exists $gene->{$unit[3]} or $unit[3]=~/sb/i);
    push @{$hash{$unit[0]}}, [$unit[3],$unit[1]];
}
close IN;
my %order;
foreach my $chr (keys %hash){
    my @sortarray=sort {$a->[1] <=> $b->[1]} @{$hash{$chr}};
    for(my $i=0;$i<@sortarray;$i++){
         $order{$sortarray[$i][0]}=[$chr,$i];
         #print "$sortarray[$i][0]\t$chr\t$i\n";
    }
}
return \%order;
}



sub parseblock
{
my ($block,$bed,$blast)=@_;
my ($refer,$refn,$genen);
my (%ref,%qry,$pair,$oryza);
my (%refgene,%qrygene);
my (%lose1,%lose2);
open OS, ">$block.os.bestdistance" or die "$!";
open OB, ">$block.ob.bestdistance" or die "$!";
open OUT, ">$block.table" or die "$!";
open OUT1, ">$block.geneinblock" or die "$!";
open OUT2, ">$block.summary";
open IN, "$block" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my $line=$_;
     if ($line=~/## View (\d+)\: pivot (\w+)$/){
        $refn=$1;
        $refer=$2; 
     }elsif($line=~/^#/){
        next;
     }elsif($line=~/^\s+(\d+)\-\s*(\d+)\:\s+?(.*)$/){
        my $genen=$2;
        my $genes=$3;
        my @unit=split("\t",$genes);
        if ($unit[0]=~/\w+/ and $unit[1]=~/\w+/){ ##syntenic gene pair
            $pair++;
            $oryza++ unless ($unit[2]=~/\w+/); ## oryza specific gene pair
            #print OUT "$refer\t$refn\t$genen\t$unit[0]\t$unit[1]\n";
            my @gene1=split(";",$unit[0]);
            my @gene2=split(";",$unit[1]);
            my $size1=@gene1;
            my $size2=@gene2;
            ##chr gene1 gene1pos tandemnumber gene2 gene2pos tandemnumber
            if ($size1 > 1 or $size2 > 1){
               print OUT "$refer\t$gene1[0]\t$bed->{$gene1[0]}\t$size1\t$gene2[0]\t$bed->{$gene2[0]}\t$size2\n";
            }
            foreach my $g1 (@gene1){
               $ref{$g1}=1;
            }
            foreach my $g2 (@gene2){
               $qry{$g2}=2;
            }
        }elsif($unit[0]=~/\w+/ and $unit[2]=~/\w+/){  ## lose in brachyantha
            my @gene1=split(";",$unit[0]);
            foreach my $g1 (@gene1){
               $refgene{$g1}=1;
               $lose2{$unit[0]}++;
            }
        }elsif($unit[1]=~/\w+/ and $unit[2]=~/\w+/){  ## lose in rice
            my @gene2=split(";",$unit[1]);
            foreach my $g2 (@gene2){
               $qrygene{$g2}=1;
               $lose1{$unit[1]}++;
            }
        }elsif($unit[0]=~/\w+/){ ## duplicate in rice
            my @gene1=split(";",$unit[0]);
            foreach my $g1 (@gene1){
              $refgene{$g1}=1;
              my $best=$blast->{$g1}->[0] ? $blast->{$g1}->[0] : "NA";
              ## distance to best hit gene, 10000 mean present on different chromosome
              #print "$g1\t$bed->{$g1}->[0]\t$best\t$bed->{$best}->[0]\n";
              my $dist=$bed->{$g1}->[0] eq $bed->{$best}->[0] ? abs ($bed->{$g1}->[1]-$bed->{$best}->[1]) : 10000;
              print OS "$g1\t$best\t$dist\n";
            }
        }elsif($unit[1]=~/\w+/){ ## duplicate in brahchyantha
            my @gene1=split(";",$unit[1]);
            foreach my $g1 (@gene1){
              $refgene{$g1}=1;
              my $best=$blast->{$g1}->[0] ? $blast->{$g1}->[0] : "NA";
              ## distance to best hit gene, 10000 mean present on different chromosome
              my $dist=$bed->{$g1}->[0] eq $bed->{$best}->[0] ? abs ($bed->{$g1}->[1]-$bed->{$best}->[1]) : 10000;
              print OB "$g1\t$best\t$dist\n";
            }
        }
     } 
}
close IN;
close OUT;
close OUT1;
close OS;
close OB;
my $refgenen=keys %ref; ## gene in pair
my $qrygenen=keys %qry; ## gene in pair
my $refgenet=$refgenen+keys %refgene; ## gene not in pair
my $qrygenet=$qrygenen+keys %qrygene; ## gene not in pair
my $lose1n=keys %lose1;
my $lose2n=keys %lose2;
#print OUT2 "Total Pair: $pair\n";
#print OUT2 "Reference Gene Number: $refgenen\n";
#print OUT2 "Target Gene Number: $qrygenen\n";
#print OUT2 "Total Ref Gene Number: $refgenet\n";
#print OUT2 "Total Qry Gene Number: $qrygenet\n";
print OUT2 "Chr\t#Pair\t#RefGeneNumber\t#TarGeneNumber\t#TotalRefGeneNumber\t#TotalTarGeneNumber\t#LoseInOs\t#LoseInOb\t#Oryza\n";
print OUT2 "$chr\t$pair\t$refgenen\t$qrygenen\t$refgenet\t$qrygenet\t$lose1n\t$lose2n\t$oryza\n";
close OUT2;
}## end of sub


