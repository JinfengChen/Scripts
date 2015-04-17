#!/usr/bin/perl

=head1 Name

rank_sum_test.pl  --  perform rank sum test for group kaks on phylogeny tree

=head1 Description

The programs invoked: R wilcox.test().
Use rank sum test to evaluate the over-representation of big dN/dS between two subtrees

=head1 Version
  
  Author: Yujie Hu, huyj@genomics.org.cn
  Modfiy: Fan Wei, fanw@genomics.org.cn
  Version: 1.1,  Date: 2008-8-2

=head1 Usage
   
   perl rank_sum_test.pl [options] <nhx_file> <kaks_file>
   --R               set the path of R software
   --minsam          set the mimimum cutoff for sample size, default=6
   --outdir          set the output directory, default "./"
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple

   perl ../bin/rank_sum_test.pl ../input/leptin.cds.fa.muscle.nj.nhx ./leptin.cds.fa.muscle.axt.kaks 

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use lib "$Bin/../lib";
use Tree::nhx;

my ($R,$Minsam,$Outdir,$Verbose,$Help);
GetOptions(
	"R:s"=>\$R,
	"minsam:i"=>\$Minsam,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if ($Help);

my $nhx_file = shift;
my $kaks_file = shift;
my $nhx_core = basename($nhx_file);
my %kaks;
my (%root,%left,%right);

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

$R ||= "/share/raid1/genome/bin/R";
$Minsam ||= 5;

##read kaks file into memory 
open (IN,"$kaks_file")||die "fail open $kaks_file";
<IN>;
while (<IN>) {
	chomp;
	my @array = split;
	if ($array[0]=~/(\S+)&(\S+)/){ 
		my $key = ($1 lt $2)?"$1&$2":"$2&$1";
		if ($array[4]=~/[-\d\.eE]+/) {
			$kaks{$key} = $array[4];
		}
	}
}
close IN;

##split the tree into 3 subtrees at each internal node
my $tree = new Tree::nhx;
$tree->parse($nhx_file,"file");
open (OUT,">$Outdir/$nhx_core.split")||die "fail open $Outdir/$nhx_core.split";
print OUT $tree->split_tree();
close OUT;

##construct relationships for nodes and genes
open (IN, "$Outdir/$nhx_core.split")||die "fail open $Outdir/$nhx_core.split\n";
$/=">";<IN>;$/="\n";
while (<IN>) {
	my @temp = split (/:/,$_);
	my $name = $temp[1];
	$name =~s/\s+$//;$name =~s/^\s+//;;
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$/="\n";

	my @array = split (/\n/,$seq);
	my @a = split (/:/,$array[0]); my $root = $a[1];
	my @b = split (/:/,$array[1]); my $left = $b[1];
	my @c = split (/:/,$array[2]); my $right = $c[1];
	$root =~s/\s+$//; $root =~s/^\s+//; 
	$left =~s/\s+$//; $left =~s/^\s+//; 
	$right =~s/\s+$//; $right =~s/^\s+//; 
	
	$root{$name}=$root;
	$left{$name}=$left;
	$right{$name}=$right;
}
close IN;

##generate the R shell
open (OUT,">$Outdir/$nhx_core.R")||die "fail open $Outdir/$nhx_core.R\n";
foreach my $n  (keys %left) {
	my @l = split (/\s+/,$left{$n}); my $l_num = @l;
	my @ri = split (/\s+/,$right{$n}); my $ri_num = @ri;
	my @ro = split (/\s+/,$root{$n}); my $ro_num = @ro;
	my $ro_l = $l_num+$ro_num;
	my $ro_ri = $ri_num+$ro_num;
	my $l_ri = $l_num+$ri_num;
	
	if($ro_l >= $Minsam && $ro_num > 0 && $l_num > 0) {
		print OUT  $n."_RootLeft_group1=c(";
		my @group1_kaks;
		for (my $i=0;$i<$ro_num ;$i++) {
			for (my $j=0;$j<$l_num ;$j++) {
				my $title = ($ro[$i] lt $l[$j])?"$ro[$i]&$l[$j]":"$l[$j]&$ro[$i]";
				push @group1_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}
		my $kaks_str = join(",",@group1_kaks);
		print OUT $kaks_str, ")\n";
		
		my @group2_kaks;
		print OUT $n."_RootLeft_group2=c(";
		for (my $i=0;$i<$ro_num-1 ;$i++) {
			for (my $j=$i+1;$j<$ro_num ;$j++) {
				my $title = ($ro[$i] lt $ro[$j])?"$ro[$i]&$ro[$j]":"$ro[$j]&$ro[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		for (my $i=0;$i<$l_num-1 ;$i++) {
			for (my $j=$i+1;$j<$l_num ;$j++) {
				my $title = ($l[$i] lt $l[$j])?"$l[$i]&$l[$j]":"$l[$j]&$l[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $kaks_str = join(",",@group2_kaks);
		print OUT $kaks_str,")\n";
		print OUT "wilcox.test(".$n."_RootLeft_group1,".$n."_RootLeft_group2,alternative ="."\"greater\")\n";
	}
	
	if($ro_ri >= $Minsam && $ro_num > 0 && $ri_num > 0) {
		print OUT $n."_RootRight_group1=c(";
		my @group1_kaks;
		for (my $i=0;$i<$ro_num ;$i++) {
			for (my $j=0;$j<$ri_num ;$j++) {
				my $title = ($ro[$i] lt $ri[$j])?"$ro[$i]&$ri[$j]":"$ri[$j]&$ro[$i]";
				push @group1_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $kaks_str = join(",",@group1_kaks);
		print OUT $kaks_str, ")\n";

		my @group2_kaks;
		print OUT $n."_RootRight_group2=c(";
		for (my $i=0;$i<$ro_num-1 ;$i++) {
			for (my $j=$i+1;$j<$ro_num ;$j++) {
				my $title = ($ro[$i] lt $ro[$j])?"$ro[$i]&$ro[$j]":"$ro[$j]&$ro[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		for (my $i=0;$i<$ri_num-1 ;$i++) {
			for (my $j=$i+1;$j<$ri_num ;$j++) {
				my $title = ($ri[$i] lt $ri[$j])?"$ri[$i]&$ri[$j]":"$ri[$j]&$ri[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $kaks_str = join(",",@group2_kaks);
		print OUT $kaks_str,")\n";
		print OUT "wilcox.test(".$n."_RootRight_group1,".$n."_RootRight_group2,alternative ="."\"greater\")\n";
	}

	if($l_ri >= $Minsam && $l_num > 0 && $ri_num > 0) {
		my @group1_kaks;
		print OUT $n."_LeftRight_group1=c(";
		for (my $i=0;$i<$ri_num ;$i++) {
			for (my $j=0;$j<$l_num ;$j++) {
				my $title = ($ri[$i] lt $l[$j])?"$ri[$i]&$l[$j]":"$l[$j]&$ri[$i]";
				push @group1_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $kaks_str = join(",",@group1_kaks);
		print OUT $kaks_str, ")\n";

		my @group2_kaks;
		print OUT $n."_LeftRight_group2=c(";
		for (my $i=0;$i<$ri_num-1 ;$i++) {
			for (my $j=$i+1;$j<$ri_num ;$j++) {
				my $title = ($ri[$i] lt $ri[$j])?"$ri[$i]&$ri[$j]":"$ri[$j]&$ri[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		for (my $i=0;$i<$l_num-1 ;$i++) {
			for (my $j=$i+1;$j<$l_num ;$j++) {
				my $title = ($l[$i] lt $l[$j])?"$l[$i]&$l[$j]":"$l[$j]&$l[$i]";
				push @group2_kaks, $kaks{$title} if ($kaks{$title} ne '');
			}
		}	
		my $kaks_str = join(",",@group2_kaks);
		print OUT $kaks_str,")\n";
		print OUT "wilcox.test(".$n."_LeftRight_group1,".$n."_LeftRight_group2,alternative ="."\"greater\")\n";
	}
}
close OUT;

##run the R shell
`$R --vanilla --slave <$Outdir/$nhx_core.R >$Outdir/$nhx_core.R.test`;

