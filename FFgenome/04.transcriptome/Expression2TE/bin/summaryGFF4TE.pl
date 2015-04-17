#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","bed:s","status:s","help");


my $help=<<USAGE;
Summary TE information from gff table. Also include distance to gene and methylation status.
--gff: TE annotation gff file
--bed: gene annotation bed file
--status: methylation status
TE261982        0.918032786885246
TE249509        1
perl $0 --gff OBa.manual.TE.gff --bed OBa.mRNA.bed --status OBa.all.manual.TE.chr.gff.Me.status
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %family=(
      "0" => "LINE",
      "1" => "SINE",
      "2" => "Copia",
      "3" => "Gypsy",
      "4" => "LTR",
      "5" => "Helitron",
      "6" => "Stowaway",
      "7" => "Tourist",
      "8" => "Mariner",
      "9" => "Harbinger",
      "10" => "hAT",
      "11" => "En-Spm",
      "12" => "MuDR",
      "13" => "DNA",
      "14" => "Unknown"
);

my $cutoff=0.1; ### methylation level cutoff
my $me;
my $unme;
my $hash;
my $refmeth=mestatus($opt{status});
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
system("$BEDtools/closestBed -a $opt{gff} -b $opt{bed} > TEsummary.closestBED");
my $refdist=closestBED("TEsummary.closestBED");
`rm TEsummary.closestBED`;

open IN, "$opt{gff}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^#/);
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my ($id,$class);
    if ($unit[8]=~/ID=(\w+?);.*?;Class=(.*?);/){
       $id=$1;
       $class=$2;
    }
    my $len =$unit[4]-$unit[3]+1;
    my $dist=$refdist->{$id};
    #my $dist=1000;
    my $status=$refmeth->{$id};
    #print "$id\t$status\t$dist\n";
    $hash=infor($class,$dist,$len,$hash);
    if ($status >= $cutoff){
        $me=infor($class,$dist,$len,$me);
    }elsif($status ne "NA"){
        $unme=infor($class,$dist,$len,$unme);
    }
}
close IN;

foreach(sort {$a <=> $b} keys %$hash){
   my $type=$family{$_};
   my $copy=$hash->{$_}->[0];
   my $len =$hash->{$_}->[1];
   my $mean=int ($len/$copy);
   my $lenkb=int ($len/1000);
   my $dist=mean($hash->{$_}->[2]);
   my $mecopy=$me->{$_}->[0];
   my $medist=mean($me->{$_}->[2]);
   my $unmecopy=$unme->{$_}->[0];
   my $unmedist=mean($unme->{$_}->[2]);
   my $methpercent=int (100*$mecopy/($mecopy+$unmecopy));
   print "$type\t$copy\t$mean\t$lenkb\t$dist\t$mecopy/$unmecopy\t$methpercent\t$medist\t$unmedist\n";
}

##########################################
sub mestatus
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}
###########################################
sub closestBED
{
#### get the distance from TE to nearest gene
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $distance;
    my $id;
    if ($unit[10] > $unit[4]){
       $distance=$unit[10]-$unit[4]+1;
    }elsif($unit[11] < $unit[3]){
       $distance=$unit[3]-$unit[11]+1;
    }else{
       $distance=0;
    }
    if ($unit[8]=~/ID=(.*?);/){
       $id=$1;
    }
    $hash{$id}=$distance;
}
close IN;
return \%hash;
}

###########################################
sub mean
{
my ($num)=@_;
my $loop=0;
my $total=0;
@$num = sort {$a <=> $b} @$num;
foreach  (@$num) {
        next if ($_ eq "NA");
        my $temp=$_;
        $total+=$temp;
        $loop++;
}
my $number=$loop;
return 0 if ($number==0);
my $mean=$total/$number;
$median=$num->[int $number/2];
return $median;
}


###########################################
sub infor
{
my ($class,$dist,$len,$hash)=@_;
    my $index;
    if ($class=~/LINE/i){
       $index=0;
    }elsif($class=~/SINE/i){
       $index=1;
    }elsif($class=~/LTR/i){
       if ($class=~/Copia/i){
          $index=2;
       }elsif($class=~/Gypsy/i){
          $index=3;
       }else{
          $index=4;
       }
    }elsif($class=~/Helitron/i){
       $index=5;
    }elsif($class=~/DNA/){
        if($class=~/Stowaway/i){
           $index=6;
        }elsif($class=~/Tourist/i){
           $index=7;
        }elsif($class=~/Mariner/i or $class=~/Tc1/i){
           $index=8;
        }elsif($class=~/Harbinger/i or $class=~/PIF/i){
           $index=9;
        }elsif($class=~/hAT/i){
           $index=10;
        }elsif($class=~/En-Spm/i){
           $index=11;
        }elsif($class=~/MuDR/i){
           $index=12;
        }else{
           $index=13;
        }
    }else{
        $index=14;
    }
    $hash->{$index}->[0]++;
    $hash->{$index}->[1]+=$len;
    push (@{$hash->{$index}->[2]},$dist);
    return $hash;
}
