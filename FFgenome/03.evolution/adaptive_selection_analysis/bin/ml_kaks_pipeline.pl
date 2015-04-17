#!/usr/bin/perl

=head1 Name

  ml_kaks_pipeline.pl  -- put pairwise kaks based on ML method onto tree branches

=head1 Description

  Inputs: 1. Alignments of MFA format; 2. NHX tree file. 
  Both can be obtained via phylogeny_pipeline.pl, the *.mfa and *.nhx file.

=head1 Version
  
  Author: Sun Ming'an <sunma@genomics.org.cn>
  Modify: fanw@genomics.org.cn
  Version: 1.3   Date: 2008-10-06
  
=usage

perl ml_kaks_pipeline.pl [options] <MFA file> <NHX file>
  --outdir   Set output directory
  --verbose  Show running progress information
  --help     Show help information

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename qw(basename);
use Data::Dumper;
use Tree;
use Tree::nhx_svg;

##Global variables
my ($Outdir, $Verbose, $Help);
my (%config, $codeml);
my ($MfaFile, $MfaCore, $PhyFile, $nhxFile, $nhxCore, $OutFile, $nhFile);
my ($KaKsLines, @W, @dN, @dS);

##get options from command line
GetOptions(
	"outdir:s"   => \$Outdir,
	"verbose"    => \$Verbose,
	"help"       => \$Help
);

if (@ARGV < 1 || $Help){
print <<"END.";
  Usage: perl $0 [options] <MFA file> <NHX file>
  -outdir   Set output directory
  -verbose  Show running progress information
  -help     Show help information
END.
exit;
}

##set out directory
$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless (-d $Outdir);

##set file names
$MfaFile = shift;
$MfaCore = basename($MfaFile);
$PhyFile = "$Outdir/$MfaCore.phy";
$nhxFile = shift;
$nhxCore = basename($nhxFile);
$nhFile  = "$Outdir/$nhxCore.nh";

##config program 
parse_config("$Bin/config.txt",\%config);
$config{$codeml} ||= 'codeml';

##change result from Mfa to phylip format, so can be used by codeml and pamp respectively
MfaFilter($MfaFile, "$Outdir/$MfaCore.filter.mfa");	
Mfa2Phy("$Outdir/$MfaCore.filter.mfa", $PhyFile);	

##change tree file from NHX to NH format
ModifyNHX($nhxFile,$nhFile);

#codeml ananylis, codeml.ctl should be put under current directory
$OutFile = "$Outdir/$MfaCore.nhx.codeml";
RunCodeml($PhyFile, $nhFile, $OutFile);
print STDERR "Codeml analysis complete.\n\n" if(defined $Verbose);

open(IN,"$OutFile")||die("can not open $OutFile\n");
while (<IN>) {
	if (/dN & dS/) {
		<IN>;<IN>;<IN>;
		$/="\n\n";
		$KaKsLines=<IN>;
		chomp $KaKsLines;
		$/="\n";
		last;
	}
}
close(IN);

my @temp=split(/\n/,$KaKsLines);
foreach  (@temp) {
	if (/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$/) {
		my $W  = $1;
		my $dN = $2;
		my $dS = $3;
		push @W,$W;
		push @dN,$dN;
		push @dS,$dS;
	}
}

my $tree1=new Tree($nhFile,"file");

## judege whether there is dis-consitent in numbers
my $num_branch=$tree1->{exter}+$tree1->{inter}-1;
if (@W != $num_branch || @dN != $num_branch || @dS != $num_branch) {
	die ("the number Ka/Ks values is not consistent with the number of branches\n");
}

##add annotation to nhx tree
##modify by fanw on 2008-9-4
$tree1->add_anno("Dn",\@dN);
$tree1->add_anno("Ds",\@dS);
$tree1->add_anno("W",\@W);

##output the tree to file
open(OUT,">$Outdir/$nhxCore.codeml.nhx")||die"Cannot open $nhxCore.codeml.nhx\n";
print OUT $tree1->output_tree();
close OUT;

##draw SVG figure
##modify by fanw on 2008-9-4
##modify by Sun Ming'an on 2008-9-22
my $nhx_svg = Tree::nhx_svg->new('show_W',2,'show_B',0,"show_ruler",1,"dist_type","dm", "width",640,"skip",20,"is_real",1);
$nhx_svg->parse("$Outdir/$nhxCore.codeml.nhx","file");
#$nhx_svg->mark_tree();

#open (OUT,">$Outdir/$nhxCore.codeml.nhx")||"fail to create $Outdir/$nhxCore.out.nhx\n";
#print OUT $nhx_svg->string_nhx_format();
#close OUT;

open OUT,">$Outdir/$nhxCore.codeml.svg" || die "fail to creat final svg figure\n";
print OUT $nhx_svg->plot;
close OUT;


############################################################################
######################## subroutines  ######################################
############################################################################

##parse the config.txt file, and check the path of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	#die "$conifg_file not exist" unless(-f $conifg_file);
	#my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
			#if (! -e $software_address){
			#	warn "Non-exist:  $software_name  $software_address\n";
			#	$error_status = 1;
			#}
		}
	}
	close IN;
	#die "\nExit due to error of software configuration\n" if($error_status);
}

###remove annotations of NHX tree, so it can be use by codeml and pamp
#######################################################
sub ModifyNHX{
	my ($in,$out) = @_;
	open(IN,"$in")||die"Cannot open $in\n";
	open(OUT,">$out")||die"Cannot open $out\n";
	while(<IN>){
		s/\[.*\]//g;
		print OUT $_;
	}
	close IN;
	close OUT;
}

### filter muscle result
################################################################################
sub MfaFilter{
	my ($MfaFile, $FilteredMfaFile) = @_;
	my (%seq, %filteredSeq, $tag, $spNum, $seqLen, $filteredSeqLen);
	$spNum=$seqLen=$filteredSeqLen=0;
	open(IN,$MfaFile)||die"Cannot open $MfaFile\n";
	while(<IN>){
		chomp;
		if(/>(\S+)/){
			$tag = $1;
			$spNum++;
		}
		else{
			$seq{$tag}.=$_;
			if($spNum==1){
				$seqLen+=length($_);
			}
		}
	}
	close IN;

	foreach my $tmp (keys %seq){
		if(length($seq{$tmp}) != $seqLen){
			die"The length of $tmp is not consistent with others\n";
		}
	}

	for(my $i=0; $i< $seqLen; $i+=3){
		my $codonStatus=1;
		my %codons = undef;
		foreach my $tmp (keys %seq){
			my $codon = substr($seq{$tmp},$i,3);
			$codons{$tmp} = $codon;
			unless($codon =~ /[ATGC]{3}/i){
				$codonStatus = 0;
				last;
			}
		}
		if($codonStatus ==1){
			$filteredSeqLen+=3;
			foreach my $tmp (keys %seq){
				$filteredSeq{$tmp} .= $codons{$tmp};
			}
		}
		elsif($codonStatus ==0){
			next;
		}
	}
	
	open(OUT,">$FilteredMfaFile")||die"Cannot open $FilteredMfaFile\n";
	foreach my $tmp (keys %filteredSeq){
		print OUT ">$tmp\n$filteredSeq{$tmp}\n";
	}
	close OUT;
	
}

##change muscle result from fasta format to phylip format
############################################################
sub Mfa2Phy {
	my ($MfaFile, $phylipFile) = @_;
	my $seqCount = 0;
	my $seq = my $seqName = "";
	open(IN, $MfaFile)||die"Couldn't open $MfaFile\n";
	while (my $line = <IN>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		if ($line =~ /^>/) {
			$seqCount++;
		}elsif ($seqCount == 1) {
			$seq .= $line;
		}
	}
	close IN;
	my $seqLen = length $seq;
	
	open(IN, $MfaFile)||die"Can't open $MfaFile\n";
	open(OUT, ">$phylipFile")||die "Can't open $phylipFile\n";
	print OUT $seqCount," ",$seqLen,"\n";
	$seqCount = 0;
	$seq = "";
	while(my $line = <IN>) {
		chomp $line;	
		next if($line =~ /^\s*$/);
	
		if($line =~ /^>(\S+)/) {
			if ($seqCount) {
				my $len = length $seq;
				if ($len == $seqLen) {
					print OUT "$seqName  $seq\n";
					$seq = $seqName = "";
				}else {
#					unlink $MfaFile;
#					unlink $phylipFile;
					die "Error: the sequence length of $seqName is not same as others.\n";
				}
			}	
			$seqName = $1;
			$seqCount++;
		}else {
			$seq .= $line;		
		}		
	}
	close IN;
	# check the length of last sequence
	my $len = length $seq;
	if ($len == $seqLen) {
		print OUT "$seqName  $seq\n";
	}else {
#	unlink $unixFile;
#		unlink $phylipFile;
		die "Error: the sequence length of $seqName is not same as others.\n";
	}	
	close IN;
	close OUT;
}

##Make codeml.ctl, then run codeml
###################################################
sub RunCodeml{
	my ($SeqFile, $TreeFile, $OutFile) = @_;
#die"There are codeml.ctl exist\n" if (-e 'codeml.ctl');
	open(OUT,">codeml.ctl")||die"Cannot creat codeml.ctl\n";
	printf OUT <<"END.";
      seqfile = $SeqFile
     treefile = $TreeFile
	  outfile = $OutFile

        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        model = 1   * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 4.54006   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)

  fix_blength = 1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0   * 0: simultaneous; 1: one branch at a time
END.
	close OUT;
	`codeml > $Outdir/codeml.log`;  # run codeml
	`rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc lnf rst rst1 rub`; # clear tmp files produced by codeml
}
