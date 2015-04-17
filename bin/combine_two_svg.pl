#!/usr/bin/perl

=head1 Name

combine two svg files

=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

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

my $svg_file1 = shift;
my $svg_file2 = shift;

my ($content1,$content2);
$content1 = read_file($svg_file1);
$content2 = read_file($svg_file2);

$content1 =~ s/<\/svg>\s*$//;
$content2 =~ s/<\?xml.+?\n<\!DOCTYPE.+?\n<svg height.+?\n//;

print $content1.$content2;


####################################################
################### Sub Routines ###################
####################################################

sub read_file {
	my $file = shift;
	my $content;
	open IN,$file || die "fail";
	while (<IN>) {
		$content .= $_;
	}
	close IN;
	return $content;
}