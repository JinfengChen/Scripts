#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"repeat:s","embl:s","title:s","help");


my $help=<<USAGE;
add gff format of repeat annotation to embl file
Run: perl $0 -repeat FF.repeatmasker.gff -embl test.embl -title FF

USAGE

if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %repeat;
open IN, "$opt{repeat}" or die "can not open my infile";
while(<IN>){
  chomp $_;
  next if ($_ eq "");
  next if ($_=~/^#/);
  my @unit=split("\t", $_);
  my $name;
  my $family;
  if ($unit[8]=~/Target=(.*)? \d+ \d+;Class=(.*);?/ or $unit[8]=~/Target=(.*)?\:\d+\.\.\d+;Class=(.*);?/){
      $name=$1;
      $family=$2;
  } 
  if($family=~/LTR/){
     $part1="FT   LTR             $unit[3]\.\.$unit[4]\n";
     $part2="FT                   /rpt_family=\"$family\"\n";
     $part3="FT                   /rpt_name=\"$name\"\n";
     $part="$part1$part2$part3";
     $repeat{$unit[3]}=$part;
  }else{
     $part1="FT   repeat_region   $unit[3]\.\.$unit[4]\n";
     $part2="FT                   /rpt_family=\"$family\"\n";
     $part3="FT                   /rpt_name=\"$name\"\n";
     $part="$part1$part2$part3";
     $repeat{$unit[3]}=$part;
  }
}
close IN;

my $genenumber=0;
my @gene;
open EMBL, "$opt{embl}" or die "can not open my infile2";
while (<EMBL>) {	
	if ($_=~/\smRNA\s/) {
           $genenumber++;
	   if (length $gene[$genenumber] > 0) {
               $gene[$genenumber].=$_;
	   }else{
	       $gene[$genenumber]=$_;
           }
	}elsif($_=~/^SQ\s/){
	   $genenumber++;
	   if (length $gene[$genenumber] > 0) {
                $gene[$genenumber].=$_;
           }else{
		   $gene[$genenumber]=$_;
           }
	}else{
           if (length $gene[$genenumber] > 0) {
               $gene[$genenumber].=$_;
           }else{
               $gene[$genenumber]=$_;
	   }
	}
}
close EMBL;
$l=@gene;
$header=shift @gene;
$tail=pop @gene;

my %cds;
foreach (@gene) {
    print "$_\n";
    if ($_=~/join\((\d+)\.\./) {
        $cds{$1}=$_;
    }else{
        print "Not pass $_\n";
    }
}

%hash=(%repeat,%cds);
my @start=sort {$a <=> $b} keys %hash;

open ME, ">$opt{title}.merge" or die "can not open my outfile";
print ME "$header";
foreach  (@start) {
	print ME "$hash{$_}";
}
print ME "$tail";
close ME;
