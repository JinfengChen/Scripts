#!/usr/bin/perl

=head1 Name

get orf for glimmer training

=head1 Description

Read in the genome or unigene file in multiple fasta format, and output 
is the orf set for build-icm program.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl glimmer_train_orf.pl <sequence.fa>
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

my $file = shift;
my $file_base = basename($file);

my $softpath = "/opt/blc/glimmer-302/bin";

my %Seq;
my $output;

Read_fasta($file,\%Seq);

foreach my $seq_id (sort keys %Seq) {
	my $seq_str = $Seq{$seq_id}{seq};
	
	my $temp_seq_file = "./$file_base.temp";
	open OUT, ">$temp_seq_file" || die "fail $temp_seq_file";
	print OUT ">$seq_id\n$seq_str";
	close OUT;
	
	`$softpath/long-orfs -n -t 1.15 $temp_seq_file $temp_seq_file.longorf`;
	`$softpath/extract -t $temp_seq_file $temp_seq_file.longorf > $temp_seq_file.train`;
	
	my %Train;
	Read_fasta("$temp_seq_file.train",\%Train);

	foreach my $orf_id (sort keys %Train) {
		my $orf_str = $Train{$orf_id}{seq};
		my $orf_head = $Train{$orf_id}{head};
		$output .= ">$seq_id\_$orf_head\n$orf_str";
	}
	`rm $temp_seq_file*`;

}

my $orf_train_file = "./$file_base.train";
open OUT, ">$orf_train_file" || die "fail $orf_train_file";
print OUT $output;
close OUT;


####################################################
################### Sub Routines ###################
####################################################



#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my $name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		##$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}
