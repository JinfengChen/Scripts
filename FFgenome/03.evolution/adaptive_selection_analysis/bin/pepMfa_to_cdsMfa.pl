#!/usr/bin/perl

=head1 Name

pepMfa_to_cdsMfa.pl  --  convert protein alignment to cds alignment

=head1 Description

generate cds multiple alignment based on protein multiple alignment

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl pepMfa_to_cdsMfa.pl [option] <pep.mfa> <cds.fa>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/pepMfa_to_cdsMfa.pl ./leptin.cds.fa.pep.fa.muscle ../input/leptin.cds.fa > ./leptin.cds.fa.muscle

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);


my $aa_align_file = shift;
my $cds_file = shift;

my %aa;

##read the protein alignment information
open IN, $aa_align_file || die "fail open $aa_align_file\n";
$/=">"; <IN>; $/="\n";
while (<IN>) {
	my $name = $1 if(/^(\S+)/);
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	$aa{$name} = $seq;
}
close IN;

##convert protein alignment to cds alignment
my $out;
open IN, $cds_file || die "fail open $cds_file\n";
$/=">"; <IN>; $/="\n";
while (<IN>) {
	$out .= ">".$_;
	my $name = $1 if(/^(\S+)/);
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$seq =~ s/\s//g;
	$seq =~ s/---//g;
	$seq = uc($seq);
	$/="\n";
	
	my $cds;
	my $prot = $aa{$name};
	my $len_prot = length($prot);
	my $j = 0;
	for (my $i=0; $i<$len_prot; $i++) {
		my $aa = substr($prot,$i,1);
		if ($aa ne '-') {
			$cds .= substr($seq,$j,3);
			$j += 3;
		}else{
			$cds .= '---';
		}
	}
	Display_seq(\$cds);
	$out .= $cds;
}
close IN;


print  $out;

	
#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################
