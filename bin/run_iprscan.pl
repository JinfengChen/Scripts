#!/usr/bin/perl

=head1 Name

run_iprscan.pl  --  the pipeline to run iprscan

=head1 Description

this program is the pipeline to run iprscan, note that this program runs 
very slowly.

For qsub_sge way, one job take up 7.6G memory, i.e, one compute-node only hold two 
iprscan jobs at the same time.

For multi_process way, only one iprscan job can run on the host compute-noed at the same time.

Use multiple -appl flags to specify multiple applications. The possible applications
are listed below:  
	blastprodom
	fprintscan
	hmmpfam
	hmmpir
	hmmpanther
	hmmtigr
	hmmsmart
	superfamily
	gene3d
	scanregexp
	profilescan
	seg
	coils


=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-5-7
  Note:

=head1 Usage
  
  perl run_iprscan.pl [options] <proteins.fa>
  --appl <str>    set applications, this option can be used multiple times
  --cuts <int>   set the number of sequences in each cutted file, default=100
  --cpu <int>	 set the cpu number to use in parallel, default=3   
  --run <str>    set the parallel type, qsub, or multi, default=qsub
  --queue <str>  set the queue
  --outdir <str>  set the result directory, default="."
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

   perl ../bin/run_iprscan.pl -cuts 10 -cpu 20 ../input/rice_prot100.fa 

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my (@Appl,$Cuts,$Cpu,$Run,$Outdir);
my ($Verbose,$Help);
my $Queue;
GetOptions(
    	"cuts:i"=>\$Cuts,
	"appl:s"=>\@Appl,
	"cpu:i"=>\$Cpu,
	"queue:s"=>\$Queue,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Cuts ||= 100;
$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $seq_file = shift;
my $seq_file_name = basename($seq_file);

my %config;
parse_config("$Bin/config.txt",\%config);

my $iprscan = "$config{iprscan} -cli -ipr -goterms -format raw ";
my $fastaDeal = $config{fastaDeal};
my $qsub_sge = $config{qsub_sge};
my $multi_process = $config{multi_process};

##add iprscan applications
foreach  (@Appl) {
	$iprscan .= " -appl $_"
}

my $iprscan_shell_file = "$Outdir/$seq_file_name.iprscan.sh";
my @subfiles;

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

`perl $fastaDeal -cuts $Cuts $seq_file -outdir $Outdir`;
@subfiles = glob("$Outdir/$seq_file_name.cut/*.*");

##creat shell file
open OUT,">$iprscan_shell_file" || die "fail $iprscan_shell_file";
foreach my $subfile (@subfiles) {
	print OUT "$iprscan -i $subfile -o $subfile.iprscan; \n";
}
close OUT;

##run the shell file
my $opt;
$opt=" -queue $Queue " if defined $Queue;
`perl $qsub_sge $opt --maxjob $Cpu  --resource vf=8.1G $iprscan_shell_file` if ($Run eq "qsub");
`perl $multi_process -cpu 1 $iprscan_shell_file` if ($Run eq "multi");

##cat together the result
`cat $Outdir/$seq_file_name.cut/*.iprscan > $Outdir/$seq_file_name.iprscan`;

##extract the IPR and GO annotations for genes
`perl $Bin/iprscan_parser.pl $Outdir/$seq_file_name.iprscan -outdir $Outdir`;


####################################################
################### Sub Routines ###################
####################################################


##parse the config.txt file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		next if(/^#/);
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


