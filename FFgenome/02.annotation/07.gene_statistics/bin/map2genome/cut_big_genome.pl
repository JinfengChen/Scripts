#!/usr/bin/perl 

=head1 Name

  cut_big_genome.pl --cut big genome sequence into some subfiles.

=head1 Version

  Author:Quanfei Huang,huangqf@genomics.org.cn
  Version:1.0
  Date:2008-10-12

=head1 Usage

  perl cut_big_genome.pl [Options] FastaFile
  --outdir   set sub files dir
  --help     print this information
  --maxjobs  set the max num of sub files

=cut

use strict;
use File::Basename qw(basename);
use Getopt::Long;

my $Outdir;
my $Help;
my $Maxjobs;
GetOptions(
	"outdir:s"=>\$Outdir,
	"help"=>\$Help,
	"maxjobs:i"=>\$Maxjobs,
);

die `pod2text $0` if (@ARGV==0 or $Help);
$Outdir ||=".";
$Maxjobs ||=5;

my $Fastafile=shift;
my $Fastafile_name=basename($Fastafile);

my %sub_seq;
my $Total_len=0;
open IN,$Fastafile or die "can not open $Fastafile:$!";
$/=">";<IN>;$/="\n";
while(<IN>) {
        chomp;
        my $head=$_;
        $/=">";
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        $seq=~s/\s//g;
        $Total_len+=length($seq);
}
close IN;

my $Sub_len=int ($Total_len/$Maxjobs);
my $Cur_len=0;
open IN,$Fastafile or die "can not open $Fastafile:$!";
$/=">";<IN>;$/="\n";
while(<IN>) {
        chomp;
        my $head=$_;
        $/=">";
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        $seq=~s/\s//g;
        $Cur_len+=length($seq);
	my $sub_id=int($Cur_len/$Sub_len);
        $sub_seq{"$sub_id"}.=">$head\n$seq\n";
}

my $i=0;
my $cut_dir="$Outdir/$Fastafile_name.cut";
`rm -r $cut_dir` if (-e $cut_dir);
mkdir($cut_dir) or die "can not mkdir $cut_dir:$!";
foreach (values %sub_seq){
        open OUT,">$cut_dir/$Fastafile_name.$i" or die"$!";
        print OUT $_;
        close OUT;
        $i++;
}


