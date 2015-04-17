#!/usr/bin/perl

=head1 Name

get gene list from hcluster result

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


my $hcluster_file = shift;


open IN, $hcluster_file || die "fail open $hcluster_file\n";

while (<IN>) {
	my (%name_seq,%pair,$output);

	chomp;
	my @t = split /\s+/;
	my @gene = split /,/,$t[-1];

	foreach my $gene_id (@gene) {
		$name_seq{$gene_id} = $t[0];
	}	
	
	foreach my $first (sort keys %name_seq) {
		foreach my $second (sort keys %name_seq) {
			next if($first eq $second || exists $pair{"$second&$first"});
			$pair{"$first&$second"} = 1;
		}
	}

	foreach (sort keys %pair) {
		if (/([^&]+)&([^&]+)/) {
			$output .= "$_\t$name_seq{$1}\n";
		}
	}
	
	print $output;

}
close IN;

