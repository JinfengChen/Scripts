#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use FindBin;

#$ ort2matchlist.pl os_sb.ort > os_sb.matchlist
#Run DAGchainer
#$ ~/DAGCHAINER/run_DAG_chainer.pl -i os_sb.matchlist -Z 12 -D 20 -g 1 -A 5 &
#
##Reverse original maitch lst
#$ cat os_sb.matchlist | perl -ne 'chomp; split; print join("\t", $_[4], $_[5], $
#_[6], $_[7], $_[0], $_[1], $_[2], $_[3], $_[8],),"\n";' > os_sb.matchlist2
#
##Analyze DAG
#$ analyze_dagchain.pl os_sb.matchlist2 os_sb.matchlist.aligncoords > os_sb.analyze_dag.out
#
##Call synteny:
#$ call_synteny.pl os_sb.analyze_dag.out fake.manual 5
#
##List blocks:
#$ define_blocks.pl os_sb.ort os_sb.matchlist.aligncoords > os_sb.blocks

my ($dist_thresh, $verbose);
GetOptions(
	"verbose"		=> \$verbose,
);

my $execdir=$FindBin::Bin . "/synteny_method/scripts";
my $dagdir=$FindBin::Bin . "/DAGCHAINER";

my $ort_file=shift;
my ($name, $dir)=fileparse($ort_file, ".ort");

my $matchlist_file="$dir$name.matchlist";
my $aligncoords_file="$matchlist_file.aligncoords";
my $blocks_file="$dir$name.blocks";

foreach my $cmd (
	"$execdir/ort2matchlist2.pl $ort_file > $matchlist_file",
	"$dagdir/run_DAG_chainer.pl -i $matchlist_file -Z 12 -D 20 -g 1 -A 5",
	"$execdir/define_blocks2.pl $ort_file $aligncoords_file > $blocks_file",
) {
	0 == system($cmd) || die "FAILED to run $cmd\n";
}

1;
