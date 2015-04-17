#!/usr/bin/perl

=head1 Name

get_branch_leaves.pl  --  get all species belong to a phylogeny branch

=head1 Description

This program is to used to get all the species(genus, or other ranks) that belong to a 
phylogeny branch, using the NCBI taxomony database.
You can use this program to retrieve leaf names at any rank(example,species,genus)
of any phylogeny branch(exaple, plant or animal).

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl get_branch_leaves.pl <branch root tax_id>
  --rank <str>   set phylogeny rank, species, genus, etc, default=species   
  --type <str>   set name type, scientific, common, all, etc, default=scientific
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

 perl ../bin/get_branch_leaves.pl 3398

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use lib "$Bin/../lib";


##get options from command line into variables and set default values
my ($Rank,$Type);
my ($Verbose,$Help);
GetOptions(
	"rank:s"=>\$Rank,
	"type:s"=>\$Type,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Rank ||= "species";
$Type ||= "scientific";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $Root = shift;

my $nodes_file = "$Bin/../dat/nodes.dmp";
my $names_file = "$Bin/../dat/names.dmp ";

my %Corr; ##store correlations of all tax_ids to form a tree
my %Rank; ##store rank of all tax_ids
my %Need; ##store all the wanted tax_ids in this branch
my %Name; ##store all the wanted names in this branch

open IN, $nodes_file || die "fail $nodes_file\n";
while (<IN>) {
        if (/^(\d+)\s+\|\s+(\d+)\s+\|\s+(.+?)\s+\|/) {
                push @{$Corr{$2}},$1;
                $Rank{$1} = $3;
                #print $1,"\t",$2,"\t",$3,"\n";
        }
}
close IN;
print STDERR "\n\nread done\n" if(defined $Verbose);

recursion($Root);
print STDERR "recursion done\n" if(defined $Verbose);



open IN, $names_file || die "fail $names_file\n";
while (<IN>) {
        if ($Type =~ /all/i) {
			if (/^(\d+)\s+\|\s+(.+?)\s+\|/) {
					if (exists $Need{$1}) {
						$Name{$2} = 1;
					}
			}
		}
		if ($Type =~ /sci/i) {
			next if(!/scientific name/);
			if (/^(\d+)\s+\|\s+(.+?)\s+\|/) {
					if (exists $Need{$1}) {
						$Name{$2} = 1;
					}
			}
		}
        
}
close IN;
print STDERR "read names done\n" if(defined $Verbose);


my $out;
foreach  (sort keys %Name) {
        $out .= $_."\n";
}
print $out;

print STDERR "task finished\n\n" if(defined $Verbose);




#####################################################
sub recursion{
        my $root = shift;
        if ($Rank{$root} eq $Rank) {
                $Need{$root} = 1;
        }
        foreach my $child (@{$Corr{$root}}) {
                recursion($child);
        }
}
