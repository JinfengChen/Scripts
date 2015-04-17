#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"cluster:s","cds:s","blastm8:s","help");


my $help=<<USAGE;
Find out these gene that 1. do not cluster with other genes; 2. unique cluster in one species; 3. shared cluster with other species.
For unique and shared cluster table, col1 is gene id and col2 is family id.
perl $0 --cluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/all.cds > log 2> log2 &
--cluster: all_vs_all.blast.m8.solar.forHC.hcluster 
--cds: all.cds
--blastm8: all_vs_all.blast.m8
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 
`rm -R ClusterGeneType`;
`mkdir ClusterGeneType`;

my $refgeneid=getfastaid($opt{cds});
my $refcluster=parsecluster($opt{cluster});
#my $refblastm8=parseblastm8($opt{blastm8});
my $refblastm8;
my ($refunique,$refshared) =parseunique($opt{cluster});

my %totalgene; 
my %nocluster;
foreach my $gene (keys %$refgeneid){
    if ($gene=~/(.*)\_(\w+)$/){
        push (@{$totalgene{$2}},$gene);
    }
    unless (exists $refcluster->{$gene}){
        if ($gene=~/(.*)\_(\w+)$/){
           push (@{$nocluster{$2}},$gene);
        }
    }
}
print "Genome\tTotalGeneNumber\tNonClusterGeneNumber\tNonBlastHitGeneNumber\tUniqueFamilyGeneNumber\tSharedFamilyGeneNumber\n";
foreach(keys %nocluster){
   my $genen=@{$nocluster{$_}};
   my $total=@{$totalgene{$_}};
   my $noblast;
   print "$_\t$total\t$genen\t";
   open OUT, ">>./ClusterGeneType/$_\.noncluster.table" or die "$!";
      foreach(@{$nocluster{$_}}){
         if (exists $refblastm8->{$_}){
            print OUT "$_\t$refblastm8->{$_}\n";
         }else{
            print OUT "$_\n";
            $noblast++; ## in noncluster gene how many gene do not have blast homolog
         }
      }
   close OUT;
   print "$noblast\t";
   open OUT, ">>./ClusterGeneType/$_\.unique.table" or die "$!";
   print OUT "$refunique->{$_}->[0]";
   close OUT;
   print "$refunique->{$_}->[1]\t";
   open OUT, ">>./ClusterGeneType/$_\.shared.table" or die "$!";
   print OUT "$refshared->{$_}->[0]";
   close OUT;
   print "$refshared->{$_}->[1]\n";
}


#############################################
sub parsecluster
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   my @gene=split(",",$unit[6]);
   foreach my $gene (@gene){
      if ($gene=~/(.*)\_(\w+)$/){
         $hash{$gene}=$2;
      }
   }
}
close IN;
return \%hash;
}

sub getfastaid
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
    $hash{$head}=1;
}
$/="\n";
return \%hash;
}


sub parseblastm8
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   my $gene=$unit[0];
   unless ($hash{$gene}){
      my $nextline=<IN>;
      chomp $nextline;
      my @array=split("\t",$nextline);
      next if ($array[0] ne $gene);
      $hash{$gene}=$nextline;
   }
}
close IN;
return \%hash;
}

sub parseunique
{
my ($file)=@_;
my %hash;
my %shared;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   my @gene=split(",",$unit[6]);
   my %family;
   foreach my $gene (@gene){
      if ($gene=~/(.*)\_(\w+)$/){
         push (@{$family{$2}},$gene);
      }
   }
   my @species=keys %family;
   my $nspecies=@species;
   my $ngene;
   my $line;
   if ($nspecies == 1){
      foreach(@{$family{$species[0]}}){
           $ngene++;
           $line.="$_\t$unit[0]\n";
      }
      $hash{$species[0]}->[0].=$line;
      $hash{$species[0]}->[1]+=$ngene;
   }else{
      foreach my $g (keys %family){
           my ($n,$l);
           foreach(@{$family{$g}}){
              $n++;
              $l.="$_\t$unit[0]\n";
           }
           $shared{$g}->[0].=$l;
           $shared{$g}->[1]+=$n;  
      }
   }
}
close IN;
return \%hash,\%shared;
}

