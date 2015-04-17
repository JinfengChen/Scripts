#!/usr/bin/perl

=head1 Name

est-map-genome.pl  --  the pipeline for EST mapping to genome

=head1 Description

This pipeline use blat as the aligning tool, use pasa as the 
alignment assembly tool. The EST mapping result as well as 
the assembled unigene mapping result will be generated, seperately.

It has been tested that, if we allow 98% identity and 98% aligning rate, then
all the 100 testing cds sequences can be mapped to the genome. So cDNA mapping of the 
same species should follow this creteria: identity >= 98% && align_rate >= 98%.
For cDNA of other species mapping, the threashold can be relaxed.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
	  Quanfei Huang, huangqf@genomics.org.cn
  Version: 2.2,  Date: 2010-3-10
 
=head1 Usage

  perl est-map-genome.pl [options] <est.fa> <genome.fa>
  --dbcut <int>    set the number to cut database file, default=3
  --tophit <num>   set the number of top hits for the result, default no limitation
  --cpu <int>	   set the cpu number to use in parallel, default=3   
  --run <str>      set the parallel type, qsub, or multi, default=qsub
  --outdir <str>   set the result directory, default="./" 
  --identity <num> set identity cutoff, default 0.95
  --alignrate <num> set alignrate cutoff, default 0.95
  --prefix <str>   set a prefix name for the gene ID in gff3
  --queque <str>   set queue when qsub
  --cluster <str>  run cluster before assembly,T/F, default: T
  --verbose        output running progress information to screen  
  --help           output help information to screen 

=head1 Exmple

  perl ../bin/est_map_genome.pl ../input/est.tissue.10000.fa ../input/kaikoall_build2.rechained.fasta.over2k &

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Cpu,$Run,$Outdir,$DBcut);
my ($Identity_cutoff,$Alignrate_cutoff,$Prefix_tag,$Best_hit);
my ($Verbose,$Help);
my $Q;
my $Cluster;
GetOptions(
	"tophit:i"=>\$Best_hit,
	"identity:s"=>\$Identity_cutoff,
	"alignrate:s"=>\$Alignrate_cutoff,
	"prefix:s"=>\$Prefix_tag,
	"cpu:i"=>\$Cpu,
	"dbcut:i"=>\$DBcut,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"queue:s"=>\$Q,
	"cluster:s"=>\$Cluster,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$DBcut ||= 3;
$Identity_cutoff ||= 0.95;
$Alignrate_cutoff ||= 0.95;
$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= "./";
$Cluster ||='T';
die `pod2text $0` if (@ARGV == 0 || $Help);

my $est_file = shift;
my $genome_file = shift;

my $est_file_basename = basename($est_file);
my $genome_file_basename = basename($genome_file);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my %config;
parse_config("$Bin/config.txt",\%config);

my $map_shell_file = "$Outdir/$est_file_basename.map.shell";
my $pasa_shell_file = "$Outdir/$est_file_basename.pasa.shell";

my @subqrs; ## sub query files
my @subdbs; ## sub database files

##cut the input est file into small files
`perl $config{"fastaDeal.pl"} -cutf $Cpu $est_file  -outdir $Outdir`;
@subqrs = glob("$Outdir/$est_file_basename.cut/*.*");
`perl $Bin/cut_big_genome.pl  -maxjobs $DBcut  -outdir $Outdir $genome_file`;
@subdbs = glob("$Outdir/$genome_file_basename.cut/*.*");

## creat map shell file
my $format_convert_program = "perl $Bin/sim4_to_gff3.pl ";
$format_convert_program .= " --tophit $Best_hit " if(defined $Best_hit);
$format_convert_program .= " --prefix $Prefix_tag " if(defined $Prefix_tag);
open OUT,">$map_shell_file" || die "fail $map_shell_file";
foreach my $subfile (@subqrs) {
	foreach my $subdb (@subdbs) {
		my $subdb_base = basename($subdb);
		print OUT "$config{blat} -out=sim4 $subdb $subfile $subfile.$subdb_base.blat.sim4;  ",
		  "$format_convert_program --identity $Identity_cutoff --alignrate $Alignrate_cutoff  $subfile.$subdb_base.blat.sim4 > $subfile.$subdb_base.blat.sim4.gff; \n";
	}
}
close OUT;


## run map shell file
my $Queue;
$Queue = " -queue $Q " if defined $Q;
if ($Run eq "qsub") {
	`perl  $config{"qsub-sge.pl"}  --resource vf=0.9G  --maxjob $Cpu  $Queue $map_shell_file`;
}
if ($Run eq "multi") {
	`perl $config{"multi-process.pl"} -cpu $Cpu $map_shell_file`;
}

##cat together the result
`cat $Outdir/$est_file_basename.cut/*.blat.sim4 > $Outdir/$est_file_basename.blat.sim4`;
`cat $Outdir/$est_file_basename.cut/*.blat.sim4.gff > $Outdir/$est_file_basename.blat.sim4.gff`;

#use pasa to do assembly
open (OUT2,">$pasa_shell_file") || die "fail $pasa_shell_file";  ## by  minjiumeng
if ( $Cluster eq 'F' ){
	print OUT2 "perl $Bin/pasa_alignment_assembly.pl $config{pasa} $Outdir/$est_file_basename.blat.sim4.gff > $Outdir/$est_file_basename.blat.sim4.gff.pasa.gff";
}else{
	print OUT2 "perl $Bin/pasa_alignment_assembly_by_cluster.pl $config{pasa} $Outdir/$est_file_basename.blat.sim4.gff > $Outdir/$est_file_basename.blat.sim4.gff.pasa.gff";
}
close OUT2;


if ($Run eq "qsub") {
	`perl  $config{"qsub-sge.pl"}  --resource vf=4.5G  --maxjob 1  $pasa_shell_file`;
}
if ($Run eq "multi") {
	`perl $config{"multi-process.pl"} -cpu 1 $pasa_shell_file`;
}

####################################################
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


