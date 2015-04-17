#!/usr/bin/perl

=head1 Name

check_data_quality.pl  --  check the cds and pep to make sure correct 

=head1 Description

this program is run before the main ortholog_paralog_family pipeline,
to make sure the input raw data are correct.

There should be two input files: cds file and pep(protein) file;
The number of cds and pep should be equal, each pair with the same ID,
The product of cds translation should be equal to the protein.
Also recommond to remove those small pep length less than a cutoff;

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-5-30
  Note:

=head1 Usage
 
  perl check_data_quality.pl <cds_file.fa> <pep_file.fa>
  --len <int> set the cutoff of protein length, default 30
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl check_data_quality.pl ATH1_5chrs.cds ATH1_5chrs.pep
  perl check_data_quality.pl -len 50 ATH1_5chrs.cds ATH1_5chrs.pep

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Verbose,$Help,$Len);
GetOptions(
	"len:i"=>\$Len,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Len ||= 30;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $cds_file = shift;
my $pep_file = shift;

my (%CDS,$cds_total,%PEP,$pep_total);

`$Bin/cds2aa.pl $cds_file > $cds_file.translation`;

$cds_total = Read_fasta("$cds_file.translation",\%CDS);
$pep_total = Read_fasta($pep_file,\%PEP);

print  "CDS number: $cds_total; PEP number: $pep_total; Number not equal!\n" if($cds_total != $pep_total && $Verbose);

foreach my $gene (sort keys %CDS) {
	my $cds_p = $CDS{$gene};
	my $pep_p = $PEP{$gene};
	$cds_p->{seq} =~ s/U/X/g;
	$pep_p->{seq} =~ s/U/X/g;


	if (abs($cds_p->{len} - $pep_p->{len}) > 1) {
		print "$gene  Length not equal!\n";
		print  "$cds_p->{len}\n$pep_p->{len}\n" if($Verbose);
	}elsif($cds_p->{seq} ne $pep_p->{seq}){
		print "$gene  Sequence not equal!\n";
		print  "$cds_p->{seq}\n$pep_p->{seq}\n" if($Verbose);
	}elsif($pep_p->{len} < $Len){
		print "$gene  Less short than 30 aa!\n" if($Len);
		print  "$pep_p->{seq}\n" if($Verbose && $Len);
	}
}

print  "Check complete!" if($Verbose);

#`rm $cds_file.translation`;

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
		$seq=~s/\s//g;
		$/="\n";
		
		if (exists $hash_p->{$name}) {
			warn "name $name is not uniq in $file";
		}

		$hash_p->{$name}{head} =  $head;
		$hash_p->{$name}{len} = length($seq);
		$hash_p->{$name}{seq} = $seq;

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}
