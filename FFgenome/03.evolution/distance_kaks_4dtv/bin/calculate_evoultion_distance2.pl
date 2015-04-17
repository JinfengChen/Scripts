#!/usr/bin/perl

=head1 Name

calculate_evolution_distance.pl  --  calculate the evolution distance for a given set of genes

=head1 Description

This program calculate values such as ks, 4dtv, etc, and draw distribution figures.

When calculate KaKs, you should always be altert about the divergence between two species.
For the YN model, it is not accurate when Ks is larger than two. At this time, you should
use other models such as LPB. 

The current 4dtv are calculated with HKY substitution models (similar to YN), we will develop
the LPB corrected programs in future, so the 4dtv could be compared to Ks.

This version reads the list file in two column format.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-12-22
  Note:

=head1 Usage
  
  perl calculate_evolution_distance.pl [options]
  --list <file>  set the ortholog gene list file
  --cds1 <file>  set the first cds file
  --cds2 <file>  set the second cds file
  --model <str>  set the KaKs model: YN,LPB, default=YN 
  --outdir <str>  set the output directory
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/calculate_evoultion_distance.pl --list ../input/NC_003155.gb.pep.NC_003888.gb.pep.reciprocal.best.list --cds1 ../input/NC_003155.gb.cds  --cds2 ../input/NC_003888.gb.cds

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../../../common_bin";
use GACP qw(parse_config);

##get options from command line into variables and set default values
my ($List_file,$CDS1_file,$CDS2_file,$Model);
my ($Outdir,$Verbose,$Help);
GetOptions(
	"list:s"=>\$List_file,
	"cds1:s"=>\$CDS1_file,
	"cds2:s"=>\$CDS2_file,
	"model:s"=>\$Model,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Outdir ||= ".";
$Model ||= "YN";
die `pod2text $0` if ($Help || !$List_file || !$CDS1_file || !$CDS2_file);

my %Ortho_list;
my %CDS1_seq;
my %CDS2_seq;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my $List_base = basename($List_file);

my $config_file = "$Bin/../../../config.txt";
my $common_bin = "$Bin/../../../common_bin";

my $muscle = parse_config($config_file,"muscle");
my $kaks_caculator = parse_config($config_file,"kaks_caculator");

read_list($List_file,\%Ortho_list);

Read_fasta($CDS1_file,\%CDS1_seq);

Read_fasta($CDS2_file,\%CDS2_seq);

##do pairwise alignment and generate axt file
`rm $Outdir/$List_base.axt`;
foreach my $gene1 (sort keys %Ortho_list) {
	my $gene2 = $Ortho_list{$gene1};
	if(!exists $CDS1_seq{$gene1} || !exists $CDS2_seq{$gene2}) {
		warn "###############CDS not found###############\n";
		next;
	}
	open OUT, ">$Outdir/$gene1.$gene2.cds" || die "fail temp";
	print OUT ">$gene1\n$CDS1_seq{$gene1}\n>$gene2\n$CDS2_seq{$gene2}\n";
	close OUT;
	`perl $common_bin/cds2aa.pl $Outdir/$gene1.$gene2.cds > $Outdir/$gene1.$gene2.pep`;
	`$muscle -in $Outdir/$gene1.$gene2.pep -out $Outdir/$gene1.$gene2.pep.mfa`;
	`perl $Bin/pepMfa_to_cdsMfa.pl $Outdir/$gene1.$gene2.pep.mfa $Outdir/$gene1.$gene2.cds > $Outdir/$gene1.$gene2.cds.mfa`;
	`perl $Bin/mfa_to_axt.pl $Outdir/$gene1.$gene2.cds.mfa > $Outdir/$gene1.$gene2.cds.mfa.axt`;
	`cat $Outdir/$gene1.$gene2.cds.mfa.axt      >> $Outdir/$List_base.axt`;
	`rm $Outdir/$gene1.$gene2*`;
}

##calculate 4dtv and kaks
`perl $Bin/calculate_4DTV_correction.pl $Outdir/$List_base.axt > $Outdir/$List_base.axt.4dtv`;
`$kaks_caculator -m $Model -i $Outdir/$List_base.axt -o $Outdir/$List_base.axt.kaks`;

`perl $Bin/draw_4dtv_kaks.pl $Outdir/$List_base.axt.4dtv`;
`perl $Bin/draw_4dtv_kaks.pl $Outdir/$List_base.axt.kaks`;

##conjoin all the cds and calculate 4dtv and kaks
`perl $Bin/conjoin_axt.pl $List_base.axt > $List_base.join.axt`;
`perl $Bin/calculate_4DTV_correction.pl $List_base.join.axt > $List_base.join.axt.4dtv`;
`$kaks_caculator -m $Model -i $Outdir/$List_base.join.axt -o $Outdir/$List_base.join.axt.kaks`;

####################################################
################### Sub Routines ###################
####################################################


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/="\n>";
		my $seq = <IN>;
		chomp $seq;
		##$seq=~s/\s//g;
		$/="\n";
		
		$hash_p->{$name} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}


##SAV1025 NC_003155       +       1292825 1293904 SCO0188 NC_003888       +       177493  178503  4e-118  66.86
sub read_list {
	my $file = shift;
	my $hash = shift;

	open IN,$file || die "fail $file";
	while (<IN>) {
		chomp;
		my @t = split /\s+/;
		my $gene1 = $t[0];
		my $gene2 = $t[1];
		$hash->{$gene1} = $gene2;
	}
	close IN;
}

