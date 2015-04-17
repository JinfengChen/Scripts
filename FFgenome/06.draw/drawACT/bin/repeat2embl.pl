#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"repeatout:s","embl:s","title:s","help");


my $help=<<USAGE;
parse repeatmasker out file and add the result to embl file
Run: perl $0 -repeat FF.repeatmasker.out -embl test.embl

USAGE

if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %repeat;
open IN, "$opt{repeatout}" or die "can not open my infile";
<IN>;
<IN>;
<IN>;
while(<IN>){
  my @unit=split(" ", $_);
  if($unit[10]=~/LTR/){
     $part1="FT   LTR             $unit[5]\.\.$unit[6]\n";
     $part2="FT                   /rpt_family=\"$unit[10]\"\n";
     $part3="FT                   /rpt_name=\"$unit[9]\"\n";
     $part="$part1$part2$part3";
     $repeat{$unit[5]}=$part;
  }else{
     $part1="FT   repeat_region   $unit[5]\.\.$unit[6]\n";
     $part2="FT                   /rpt_family=\"$unit[10]\"\n";
     $part3="FT                   /rpt_name=\"$unit[9]\"\n";
     $part="$part1$part2$part3";
     $repeat{$unit[5]}=$part;
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
