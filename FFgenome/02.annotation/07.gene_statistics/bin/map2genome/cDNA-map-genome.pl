#!/usr/bin/perl

=head1 Name

cDNA-map-genome.pl  --  the pipeline for cDNA mapping to genome

=head1 Description

This pipeline use blat with -fine and -out=sim4 options for the 
aligning software, hits that pass the cutoff threshold will be
output to the result file in gff3 format. 

For cDNA of the same species, we require a cutoff: 
identity >= 95% && alignrate >= 90%; For cDNA of other species,
the cutoff can be relaxed based on the distance between the 
two species.

This pipeline will integrate function annotation, infer cds and 
check cds model, automatically. 

At last it will check the redandence and get a non-redundant cDNA
defined gene set.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.1,  Date: 2008-10-30
 
=head1 Usage
  
  perl cDNA-map-genome.pl [options] <cDNA.fa> <genome.fa>
  --dbcut <int>    set the number to cut database file, default=3
  --tophit <num>   set the number of top hits for the result, default no limiatation
  --identity <num>  set identity cutoff, default 0.95
  --alignrate <num>  set alignrate cutoff, default 0.95
  --cpu <int>	 set the cpu number to use in parallel, default=3   
  --run <str>    set the parallel type, qsub, or multi, default=qsub
  --outdir <str>  set the result directory, default="./"
  --prefix <str>    set a prefix name for the gene ID in gff3
  --queue <str>     set the queue to run the jobs if run='qsub',default=all.q
  --verbose         output running progress information to screen  
  --help            output help information to screen 

=head1 Exmple

 perl ../bin/cDNA-map-genome.pl ../input/Bombyx_mori.mrna.100  ../input/test_chr_123.seq &

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Cpu,$Run,$Outdir,$DBcut);
my ($Identity_cutoff,$Alignrate_cutoff,$Prefix_tag,$Best_hit);
my ($Verbose,$Help);
my $Queue;
GetOptions(
	"cpu:i"=>\$Cpu,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"tophit:i"=>\$Best_hit,
	"identity:s"=>\$Identity_cutoff,
	"alignrate:s"=>\$Alignrate_cutoff,
	"dbcut:i"=>\$DBcut,
	"prefix:s"=>\$Prefix_tag,
	#"queue:s"=>\$Queue,#add by huangqf
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$DBcut ||= 3;
$Identity_cutoff ||= 0.95;
$Alignrate_cutoff ||= 0.90;
$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= ".";
$Queue ||="all.q";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $cDNA_file = @ARGV[0];
my $genome_file = @ARGV[1];

my $cDNA_file_basename = basename($cDNA_file);
my $genome_file_basename = basename($genome_file);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my %config;
parse_config("$Bin/config.txt",\%config);

my $map_shell_file = "$Outdir/$cDNA_file_basename.map.shell";

my @subqrs; ## sub query files
my @subdbs; ## sub database files

##cut the input est file into small files
`perl $config{"fastaDeal.pl"} -cutf $Cpu $cDNA_file -outdir $Outdir`;
@subqrs = glob("$Outdir/$cDNA_file_basename.cut/*.*");
#`perl $config{"fastaDeal.pl"} -cutf $DBcut $genome_file  -outdir $Outdir`;
`perl $Bin/cut_big_genome.pl -maxjobs $DBcut -outdir $Outdir $genome_file`;
@subdbs = glob("$Outdir/$genome_file_basename.cut/*.*");

## creat map shell file
open OUT,">$map_shell_file" || die "fail $map_shell_file";
foreach my $subfile (@subqrs) {
	foreach my $subdb (@subdbs) {
		my $subdb_base = basename($subdb);
		print OUT "$config{blat} -fine -out=sim4 $subdb $subfile $subfile.$subdb_base.blat.sim4;\n",
	}
}
close OUT;

## run map shell file
if ($Run eq "qsub") {#do not use --reqsub,huangqf
	#`perl $config{"qsub-sge.pl"} --queue $Queue  --maxjob $Cpu $map_shell_file`;
	`perl $config{"qsub-sge.pl"} --resource vf=0.5G $map_shell_file`;
}
if ($Run eq "multi") {
	`perl $config{"multi-process.pl"} -cpu $Cpu $map_shell_file`;
}


##cat together the result
`cat $Outdir/$cDNA_file_basename.cut/*.blat.sim4 > $Outdir/$cDNA_file_basename.blat.sim4`;

##convert file format
my $format_convert_program = "perl $Bin/sim4_to_gff3.pl ";
$format_convert_program .= " --tophit $Best_hit " if(defined $Best_hit);
$format_convert_program .= " --prefix $Prefix_tag " if(defined $Prefix_tag);
`$format_convert_program  --identity $Identity_cutoff --alignrate $Alignrate_cutoff $Outdir/$cDNA_file_basename.blat.sim4 $cDNA_file > $Outdir/$cDNA_file_basename.blat.sim4.gff`;


###################################################
################### Sub Routines ###################
####################################################

##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;
			if (! -e $software_address){
				warn "Non-exist:  $software_name  $software_address\n"; 
				$error_status = 1;
			}
		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}



__END__


input mRNA file format:

>AB001052  624  76..450  [Bombyx mori]  Bombyx mori mRNA for histone H2A-like protein, complete cds.
tacggacgttcgcgagacacgcgtgcgttcgagtgctttcgtgtgttata
tcgttaacttttttaaacttcaaacatgtccggtcgcggaaaaggcggaa
aagttaagggcaaggtcaagtcccgttcgaaccgtgccggtcttcagttt
ccggtcggtcgtatacacagattgttgcgcaacggaaattacgctgaacg
cgttggtgccggtgcaccggtttacctggccgccgtcatggaatacttgg
ccgctgaagttttggaattggccggtaacgcagcaagagacaacaagaag
actagaattattcctagacatcttcaactcgccataaggaacgacgagga
actgaacaaactcccttccggtgtgacaatcgctcaaggcggagttttac
caaacattcaagcggtactactcccgaagaagaccgagaagaaagcttaa
aaaacgcttcaaactcgctcgcaaggacaacaacaacaacacaacatgtc
gtcgataatatttatatatgtataataaataataataatcgacgaacata
ttttgttcgttgttgttattattattattatattgtattgttgacaaaaa
tcaaaggcccttttcagggccgct
>AB001053  344  64..318  [Bombyx mori]  Bombyx mori mRNA for VAP-peptide, complete cds.
aaaaaacacagcacttagctcatcggcagaacacatctagtttgttattt
gaaagaccgcaaaatgttcaagttgacagtaattttcgctattatcgctg
tggcccaagcgggcgtcatagccccagtggtgcctgtagcacaccccgtc
gtggctcacacggccgtggtccacccggtcccactagtgcgcgctgccca
cgtggttcacaccgccccagtggttgccgctgccccagtggtcgccgctc
cggtggtggctgcggctcctatcgtgccgatagttaaacatgcgccaatt
atcgctgtccatcattaattgtagaaataaataaatatattttt


