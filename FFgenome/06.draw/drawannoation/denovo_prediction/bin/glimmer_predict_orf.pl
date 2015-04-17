#!/usr/bin/perl

=head1 Name

get orf for glimmer training

=head1 Description

Read in the train_cds file, build the icm models, predict genes, and convert the file formats. 

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl glimmer_predict_orf.pl <train_cds.fa> <genome_sequence.fa>
  --prefix <str>    set the gene prefix tag
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
my ($Verbose,$Help,$Prefix);
GetOptions(
	"prefix:s"=>\$Prefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $train_cds_file = shift;
my $genome_sequence_file = shift;
my $train_cds_base = basename($train_cds_file);
my $genome_sequence_base = basename($genome_sequence_file);

my $softpath = "/opt/blc/glimmer-302/bin";

`$softpath/build-icm -F -r $train_cds_base.icm <$train_cds_file`;
`$softpath/glimmer3 -Aatg -o0 -g110 -z1 -t30  -l $genome_sequence_file $train_cds_base.icm $genome_sequence_base`;
`mv $genome_sequence_base.predict $genome_sequence_base.glimmer`;
`perl $Bin/predict_convert.pl --checkcds --predict glimmer --minicds 150 --filter_ns --final $Prefix --log --verbose $genome_sequence_base.glimmer $genome_sequence_file`;

print "perl $Bin/predict_convert.pl --checkcds --predict glimmer --minicds 150 --filter_ns --final $Prefix --log --verbose $genome_sequence_base.glimmer $genome_sequence_file\n";

####################################################
################### Sub Routines ###################
####################################################

