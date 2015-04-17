#!/usr/bin/perl

=head1 Name

conjoin_axt.pl  --  join the axt file together 

=head1 Description

generate a sinlge CDS alingment in axt format

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl conjoin_axt.pl  <file.axt>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



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


####################################################
################### Sub Routines ###################
####################################################

my $axt_file = shift;

my $cds1;
my $cds2;

$/ = "\n\n";
open IN,$axt_file || die "fail $axt_file";
while (<IN>) {
	my @t = split /\n/;
	$cds1 .= $t[1];
	$cds2 .= $t[2];
}
close IN;

print "species1&species2\n$cds1\n$cds2\n\n";