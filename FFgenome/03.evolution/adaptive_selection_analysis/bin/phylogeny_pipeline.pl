#!/usr/bin/perl

=head1 Name

phylogeny_pipeline.pl  --  build phylogeny tree for genes in a family

=head1 Description

The programs invoked: cds2aa.pl , muscle , treebest, and SVG & Tree package.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-8-1

=head1 Usage
  
   perl phylogeny_pipeline.pl [option] <cds.fa>
   --equal_branch    use equal branch length
   --outdir          set the output directory, default "./"
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple

  perl ../bin/phylogeny_pipeline.pl ../input/leptin.cds.fa 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../lib";
use SVG;
use Tree::nhx_svg;

##global variables
my ($Outdir,$Equal_branch,$Verbose,$Help,$Cds_file,$Cds_file_core);
my (%config,$muscle,$treebest,$kaks_calculator,$R);

##get options from command line
GetOptions(
	"equal_branch"=>\$Equal_branch,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV < 1 || $Help);

$Cds_file = shift;
$Cds_file_core = basename($Cds_file);

$Outdir ||= ".";
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

parse_config("$Bin/config.txt",\%config);
$muscle = $config{muscle}; 
$treebest = $config{treebest};

##check and convert cds to protein
`$Bin/cds2aa.pl -check $Cds_file > $Outdir/$Cds_file_core.check`;		#check
check_cds("$Outdir/$Cds_file_core.check");
`$Bin/cds2aa.pl  $Cds_file > $Outdir/$Cds_file_core.pep.fa`;			#convert
print STDERR "cds convert to protein complete\n" if(defined $Verbose);

##protein multiple alignment and convert to cds multiple alignment
`$muscle -in $Outdir/$Cds_file_core.pep.fa -out $Outdir/$Cds_file_core.pep.fa.muscle 2> $Outdir/$Cds_file_core.pep.fa.muscle.log `;
`perl $Bin/pepMfa_to_cdsMfa.pl $Outdir/$Cds_file_core.pep.fa.muscle $Cds_file  > $Outdir/$Cds_file_core.muscle`;
print STDERR "muscle for pep and cds complete\n" if(defined $Verbose);

##build tree, the best method needs species tree and gene_id must be in treefam format
##so it is easy to use the nj method instead of the best method, we choose dm dn-ds merge (tree merge)
`$treebest nj  -t dm  $Outdir/$Cds_file_core.muscle > $Outdir/$Cds_file_core.muscle.nj.nhx  2> $Outdir/$Cds_file_core.muscle.nj.log`; 
print STDERR "build tree complete\n" if(defined $Verbose);


##draw svg figure
my $is_real = (defined $Equal_branch) ? 0 : 1;
my $nhx_svg = Tree::nhx_svg->new('show_W',2,'show_B',1,"show_ruler",1,"dist_type","dm", "skip",20, "is_real",$is_real);
$nhx_svg->parse("$Outdir/$Cds_file_core.muscle.nj.nhx","file");
open OUT,">$Outdir/$Cds_file_core.muscle.nj.nhx.svg" || die "fail creat svg";
print OUT $nhx_svg->plot;
close OUT;
print STDERR "draw svg figure complete\n" if(defined $Verbose);

`/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx -t pdf $Outdir/$Cds_file_core.muscle.nj.nhx.svg`;

####################################################
################### Sub Routines ###################
####################################################


##check the input cds file
####################################################
sub check_cds{
	my $file = shift;
	open IN, $file || die "fail";
	<IN>; ##read title line
	while (<IN>) {
		my $triple = (split /\t/)[4];
		die "the cds length is not triple" unless $triple;
	}
	close IN;
}


##parse the config.txt file, and check the path of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
#	die "$conifg_file not exist" unless(-f $conifg_file);

#	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
#			if (! -e $software_address){
#				warn "Non-exist:  $software_name  $software_address\n"; 
#				$error_status = 1;
#			}
		}
	}
	close IN;
#	die "\nExit due to error of software configuration\n" if($error_status);
}
