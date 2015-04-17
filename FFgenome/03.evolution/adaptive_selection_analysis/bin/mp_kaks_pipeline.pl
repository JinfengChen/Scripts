#!/usr/bin/perl

=head1 Name

  mp_kaks_pipeline.pl -- put pairwise kaks based on MP method onto tree branches

=head1 description

  Inputs: 1. multiple alignments of MFA format; 2. NHX tree file
  Both can be obtained via phylogeny_pipeline.pl, the *.nhx and *.mfa file.

=Version
  Author: Sun Ming'an <sunma@genomics.org.cn>
  Version: 1.3  Date: 2008.10.06

=head1 Usage

  perl mp_kaks_pipeline.pl [options] <MFA file> <NHX tree file>
  --outdir <str>  Set the output directory
  --verbose       output running progress information to screen
  --help          output help information

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;
use lib "$Bin/../lib";
use SVG;
use Tree::nhx_svg;

##Global variablest
my ($Outdir, $Verbose, $Help);
my (%config, $pamp, $kaks_calculator);
my ($MfaFile, $MfaCore, $nhxFile, $nhxCore, $PhyFile, $pampFile, $nhFile);
my (%ancseq, %extseq, $seqs, %info, %rootinfo, $treeLine, $kaksFile);

##get options from command line
GetOptions(
	"outdir:s" => \$Outdir,
	"verbose"  => \$Verbose,
	"help"     => \$Help
);

if (@ARGV < 1 || $Help){
print <<"END.";
  Usage: perl $0 [options] <MFA file> <NHX file>
  -outdir   outdir
  -verbose  verbose
  -help     help
END.
exit;
}

##set outdir
$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

##set file names
$MfaFile = shift;
$MfaCore = basename($MfaFile);
$PhyFile = "$Outdir/$MfaCore.phy";
$nhxFile = shift;
$nhxCore = basename($nhxFile);
$nhFile  = "$Outdir/$nhxCore.nh";


##set program path
parse_config("$Bin/config.txt",\%config);
$pamp = $config{pamp};
$kaks_calculator = $config{kaks_caculator};

##change result from Mfa to phylip phylip, so can be used by codeml and pamp respectively
MfaFilter($MfaFile, "$Outdir/$MfaCore.filter.mfa");	
Mfa2Phy("$Outdir/$MfaCore.filter.mfa", $PhyFile);				

#modify the tree file, so it can be used by pamp
ModifyNHX($nhxFile, $nhFile);
open(IN,$nhFile)||die;
my @arr = <IN>;
close IN;
my $spNum = 0;
foreach (@arr){
	if(/[^\)]\:\d/){
		$spNum++;
	}
}
open(OUT,">$Outdir/$nhxCore")||die"Cannot open $Outdir/$nhxCore\n";
print OUT "$spNum 1\n", @arr;
close OUT;

##run pamp
RunPamp($PhyFile,"$Outdir/$nhxCore","$Outdir/$MfaCore.pamp",0); 
#`mv pamp.ctl $Outdir`;	
ParsePamp("$Outdir/$MfaCore.pamp", "$Outdir/$MfaCore.axt");	

##calculate pairwise KaKs
`$kaks_calculator -i $Outdir/$MfaCore.axt -o $Outdir/$MfaCore.axt.YN -m YN > $Outdir/KaKs_Calculator.log`;
$kaksFile = "$Outdir/$MfaCore.axt.YN";

##parse KaKs result
parseKaKs($kaksFile,$nhFile,"$Outdir/$nhxCore.pamp.nhx");

##draw svg figure
my $nhx_svg = Tree::nhx_svg->new('show_W',2,'show_B',0,"show_ruler",1,"dist_type","dm", "width",640,"skip",20,"is_real",1);
$nhx_svg->parse("$Outdir/$nhxCore.pamp.nhx","file");
#$nhx_svg->mark_tree();

#open (OUT,">$Outdir/$nhxCore.pamp.nhx")||"fail to create $Outdir/$nhxCore.pamp.nhx\n";
#print OUT $nhx_svg->string_nhx_format();
#close OUT;

open OUT,">$Outdir/$nhxCore.pamp.svg" || die "fail to creat final svg figure\n";
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
		die "Error: the sequence length of $seqName is not same as others.\n";
	}	
	close IN;
	close OUT;
}

## Make pamp.ctl file, then run pamp
#########################################################
sub RunPamp{
	my ($SeqFile, $TreeFile, $outFile, $SeqType) = @_;
	open(CTL,">pamp.ctl")||die"Cannot write to pamp.ctl\n";
	print CTL<<" END.";
 seqfile  = $SeqFile
 treefile = $TreeFile
 outfile  = $outFile

 seqtype  = $SeqType
 ncatG    = 8
 nhomo    = 0
 END.
close CTL;
`pamp > $Outdir/pamp.log`;
}

## parse pamp result
##########################################################
sub ParsePamp{
	my ($PampFile, $axtFile) = @_;
	my ($spLines, %spTag);
	my $spNum=1;
	open(IN,"$PampFile")||die "can not open $PampFile\n";
    while(<IN>){
		if(/Frequencies../){		#for species list
			<IN>;
			$/="\n\n";
			$spLines = <IN>;
			chomp $spLines;
		}
		$/="\n";
		
		if(/\(1\) Branch lengths/){	#for ancient branches
			my $tmp = <IN>;
			while($tmp =~ /(\d+)\.\.(\d+)/g){
				if (exists $ancseq{$1}{'left'}){
						$ancseq{$1}{'right'}=$2;
						}
				else{
					$ancseq{$1}{'left'}=$2;
				}
			}
		}
			
		if(/list of extant/){		#for extant and recontrusted sequences
			while( <IN>){
				if(/node #(\d+)\s+([ATGC\s]+)/){
					my $sp = $1;
					my $seq = $2; $seq =~ s/\s//g;
					$ancseq{$sp}{'seq'}=$seq;
					}
				elsif(/^(\S+)\s+([ATGC\s]+)/){
					my $sp = $1;
					my $seq = $2; $seq =~ s/\s//g;
					$extseq{$sp} = $seq;
					}
				}
			}
	}
	close IN;

	my @lines = split("\n", $spLines);
	foreach (@lines){						#get species list
		if(/^(\S+)\s+\d/){
			$spTag{$spNum}=$1;
			$spNum++;
		}
	}
	
	open(OUT,">$axtFile")||die"Cannot open $axtFile";
	foreach my $key ( keys %ancseq){		#make AXT file			
		my ($left, $right, $leftseq, $rightseq);
		$left = $ancseq{$key}{'left'};
		$right = $ancseq{$key}{'right'};
		if(exists $extseq{$spTag{$left}}){
			$leftseq=$extseq{$spTag{$left}};
		}
		else{
			$leftseq=$ancseq{$left}{'seq'};
		}
		if(exists $extseq{$spTag{$right}}){
			$rightseq=$extseq{$spTag{$right}};
		}
		else{
			$rightseq=$ancseq{$right}{'seq'};
		}
		print OUT "$key-$left-$spTag{$left}\n$ancseq{$key}{'seq'}\n$leftseq\n\n";
		print OUT "$key-$right-$spTag{$right}\n$ancseq{$key}{'seq'}\n$rightseq\n\n";
	}
	close OUT;
}

##parse result of KaKs_Calculator
#############################################
sub parseKaKs{
	my ($kaksResult, $nhTree, $outTree) = @_;
	my $nodeKaks = 0;
	my $nodeTree = 1;

##Extract Dn, Ds, W from result of KaKs_Calculator
	open(IN,$kaksResult)||die"can't open $kaksResult\n";	
	while(<IN>){
		chomp;
		my ($root, $sp);
		my @line = split;
		next if /^Sequence/;
		if($line[0] =~ /^(\d+)-\d+-(\S+)$/){
			$root = $1;
			$sp = $2;
		}
		elsif($line[0] =~ /^(\d+)-(\d+)-$/){
			$root = $1;
			$sp = $2;
		}
		$info{$sp}{'root'} = $root;
		$info{$sp}{'Dn'} = ($line[2] eq 'NA') ? 'NA':sprintf("%.6f",$line[2]);
		$info{$sp}{'Ds'} = ($line[3] eq 'NA') ? 'NA':sprintf("%.6f", $line[3]);
		$info{$sp}{'W'} =  ($line[4] eq 'NA') ? 'NA':sprintf("%.6f",$line[4]);
		if($info{$sp}{'Dn'} ne 'NA' && $info{$sp}{'Ds'} eq 'NA')
		{
			$info{$sp}{'W'} = 'infinite';
			}
		elsif($info{$sp}{'Dn'} eq 'NA' && $info{$sp}{'Ds'} ne 'NA')
		{
			$info{$sp}{'W'} = '0.000000';
			}
		else{;}
			
	}
	close IN;

##Extracting tree info from NHX tree file	
	open(IN,$nhTree)||die "Cannot open $nhTree\n";
	while(<IN>){
		$treeLine.=$_;
	}
	close IN;
	if($treeLine =~ /^[^\(]*(\(.*\))[^\)]*$/s){		#refine the tree
			$treeLine = $1;
			}
	
##judge if the node number of NHX tree file is consistant with kaks result file
	while($treeLine =~ /\(/g){
		$nodeTree++;
		}
			
	foreach my $key (keys %info){
			if($key =~ /\D+/ ){
				$nodeKaks++;
				}
			}
	unless($nodeTree == $nodeKaks){
		print "$kaksResult and $nhTree do not matching\n";
		exit;
	}
	
	nh2kaks1(\$treeLine);
	nh2kaks2(\$treeLine);
	$treeLine.=';';
	open(OUT, ">$outTree")||die"cannot open $outTree\n";
   	print OUT $treeLine;
   	close OUT;
}
			
## leaf => node => root
###############################################################
sub nh2kaks1{
	my $pTree =  shift;
	my ($seq1,$seq1len,$seq2,$seq2len,$rootnum, $tag1, $tag2);
	if ($$pTree =~ /\((\w+:.*),\n(\w+:.*)\n\)/ ){	
		$tag1 = $1;
		$tag2 = $2;
		($seq1,$seq1len) = split (/:/,$1);
		($seq2,$seq2len) = split (/:/,$2); 
		if ($info{$seq1}{'root'} == $info{$seq2}{'root'}){
			$rootnum=$info{$seq1}{'root'};
		$rootinfo{$rootnum}= "($1\[\&\&NHX:Dn=$info{$seq1}{'Dn'}:Ds=$info{$seq1}{'Ds'}:W=$info{$seq1}{'W'}\],\n"."$2\[\&\&NHX:Dn=$info{$seq2}{'Dn'}:Ds=$info{$seq2}{'Ds'}:W=$info{$seq2}{'W'}\]\n)";
		$$pTree =~ s/\($tag1,\n$tag2\n\)/$rootnum/;
		}
		else{
			print "The root of $seq1 and $seq2 are different!!\n";
			return;
			}
	}
	else{
			return;
	}
	nh2kaks1($pTree);
	return;
}

## root => node => leaf
###################################################################################
sub nh2kaks2 {
	my $p2Tree = shift;
	my $tmp;
	if ($$p2Tree =~ /^(\d+)$/){
		$tmp=$rootinfo{$1};
		$$p2Tree =~ s/^$1/$tmp/ ;
	}
	elsif($$p2Tree =~ /\,\n(\d+):/){
		$tmp=$rootinfo{$1};
		$$p2Tree =~ s/\,\n$1:/\,\n$tmp:/ ;
	}
	elsif($$p2Tree =~ /\((\d+):/){
		$tmp = $rootinfo{$1};
		$$p2Tree =~ s/\($1:/\($tmp:/ ;
	}
	else{
		return;
	}
	nh2kaks2($p2Tree);
	return;
}
