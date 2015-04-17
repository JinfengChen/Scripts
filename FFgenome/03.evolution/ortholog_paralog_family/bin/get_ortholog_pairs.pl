#!/usr/bin/perl

=head1 Name

get ortholog pairs for treebest .ortho result

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

my $ortho_file = shift;

my $is_tag = 0;
open IN,$ortho_file || die "fail";
while (<IN>) {
	$is_tag = 0 if(/^\@end full_ortholog/);
	if ($is_tag) {
		chomp;
		my @t = split /\s+/;
		my @ary = ($t[0],$t[1]);
		@ary = sort @ary;

		print "$ary[0]\&$ary[1]\t$t[2]\t$t[3]\t$t[4]\n" if($t[4] == 1);

	}
	$is_tag = 1 if(/^\@begin full_ortholog/);
}
close IN;

