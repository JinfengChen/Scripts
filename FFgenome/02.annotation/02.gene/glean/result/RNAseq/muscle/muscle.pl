#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;
use File::Basename qw(basename dirname);

GetOptions(


);

die "Usage:$0 <option> <ortholog.list> <database.fa> <query.fa>\n" if( @ARGV < 2);

my $ortholog=shift;
my $database=shift;
my $query=shift;

my $ortholog_basename=basename($ortholog);
my $database_basename=basename($database);
my $query_basename=basename($query);

############################################################
############################################################
my %database_seq;
open (IN1,$database) || die "$!";
$/=">";<IN1>;$/="\n";
while(<IN1>){
	my $seq_name=$1 if($_=~/^(\S+)/);
	$/=">";
	my $sequence=<IN1>;
	chomp($sequence);
	$sequence=~s/\s+//g;
	$sequence=~tr/atcgu/ATCGU/;
	$database_seq{$seq_name}=$sequence;
	$/="\n";
}
close IN1;

my %query_seq;
open (IN2,$query) || die "$!";
$/=">";<IN2>;$/="\n";
while(<IN2>){
	my $seq_name=$1 if($_=~/^(\S+)/);
	$/=">";
	my $sequence=<IN2>;
	chomp($sequence);
	$sequence=~s/\s+//g;
	$sequence=~tr/atcgu/ATCGU/;
	$query_seq{$seq_name}=$sequence;
	$/="\n";
}
close IN2;

my $out_dir = "./muscle_result";
`rm -r $out_dir` if(-e "$out_dir");
mkpath("$out_dir");

my $number=1;
open (IN3,$ortholog) || die "$!";
while(<IN3>){
	my $muscle_out="$out_dir/muscle$number.fa";
	open (OUT1,">$muscle_out") || die "can't open: $!";
	chomp;
	my @c=split/\t/,$_;
	my $gene_ID_1=$c[0];
	my $gene_ID_2=$c[5];
	if(exists $database_seq{$gene_ID_1}){
		print OUT1 ">$gene_ID_1\n$database_seq{$gene_ID_1}\n>$gene_ID_2\n$query_seq{$gene_ID_2}\n";
	}elsif(exists $database_seq{$gene_ID_2}){
		print OUT1 ">$gene_ID_2\n$database_seq{$gene_ID_2}\n>$gene_ID_1\n$query_seq{$gene_ID_1}\n";
	}
	$number++;
}
close OUT1;

############################################################################
############################################################################
my $muscle_sh="./muscle.sh";
my @subfile=glob("$out_dir/*");
open (OUT2,">$muscle_sh") || die "can't open:$!";
foreach (@subfile) {
	print OUT2 "/nas/GAG_02/minjiumeng/GACP-8.0/software/muscle/muscle -in $_ -out $_.muscle.out\n";
}
close OUT2;

`perl /nas/GAG_02/minjiumeng/software2/qsub-sge.pl  --maxjob 50 -resource vf=0.9G $muscle_sh`;

	



