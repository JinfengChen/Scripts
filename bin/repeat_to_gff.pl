#!/usr/bin/perl

=head1 Name

repeat_to_gff.pl  --  convert repeat raw formats to gff format

=head1 Description

This program can read repeatmasker .out file, repeatproteinmask .annot file, or trf .dat file,
and convert it to gff format.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  --prefix <str>  set a prefix before repeat element ID
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

 perl ../bin/repeat_to_gff.pl ./rice.frag1M.fa.trf.dat
 perl ../bin/repeat_to_gff.pl ./rice.frag1M.fa.RepeatMasker.out
 perl ../bin/repeat_to_gff.pl ./rice.frag1M.fa.Proteinmask.annot

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

my $repeat_file = shift;

dat_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.dat$/);
out_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.out$/);
annot_to_gff3($repeat_file,$Prefix) if($repeat_file =~ /\.annot$/);

####################################################
################### Sub Routines ###################
####################################################


##facilitate to creat marks
####################################################
sub create_marker {
	my $number = shift || 100000;
	$number =~ s/\d/0/g;
	$number++;
	return $number;
}


##Start	End	PeriodSize 	CopyNumber	ConsensusSize	PercentMatches	PercentIndels	Score	A	C	G	T	Entropy(0-2)	consensus	repeatSequences
##19670039 19670073 4     8.8         4              93                0              61    22  0   28 48     1.51           TGTA       TGTATGTATGTATGTATGTATGTAGGTATGTATGT
####################################################
sub dat_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;
	
	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);
	my $chr;

	open IN,$file || die "fail open $file";
	while (<IN>) {
		$chr = $1 if(/^Sequence:\s+(\S+)/);
		my @t = split /\s+/;
		next if(@t != 15);
		my $start = $t[0];
		my $end = $t[1];
		my $id = $pre_tag."TR".$mark;
		
		my $score = $t[7];
		my $strand = "+";
		$output .= "$chr\tTRF\tTandemRepeat\t$start\t$end\t$score\t$strand\t.\tID=$id;PeriodSize=$t[2];CopyNumber=$t[3];PercentMatches=$t[5];PercentIndels=$t[6];Consensus=$t[13];\n";
		$mark++;
	}
	close IN;

	open OUT,">$file.gff" || die "fail creat $file";
	print OUT "##gff-version 3\n$output";
	close OUT;	

}



##  SW   perc perc perc  query     position in query              matching       repeat       position in repeat
##score   div. del. ins.  sequence  begin    end          (left)   repeat         class/family begin   end   (left)  ID
#245   35.2  2.5  0.6  chr1       1001400  1001556 (18753082) + TS2            SINE            80   239  (416)  12  
####################################################
sub out_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;
	
	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);

	open IN,$file || die "fail open $file";
	while (<IN>) {
		s/^\s+//;
		my @t = split /\s+/;
		next if($t[0] =~ /\D/ || !$t[0]);
		my $start = $t[5];
		my $end = $t[6];
		my $id = $pre_tag."TE".$mark;
		my $chr = $t[4];
		my $score = $t[0];
		my $strand = ($t[8] eq '+') ? "+" : "-";
		my $target = $t[9];
		my $class = $t[10];
		my @ary;
		push @ary,$t[11] if($t[11] !~ /[\(\)]/);
		push @ary,$t[12] if($t[12] !~ /[\(\)]/);
		push @ary,$t[13] if($t[13] !~ /[\(\)]/);
		@ary = sort {$a<=>$b} @ary;
		my ($target_start,$target_end) = ($ary[0],$ary[1]);

		$output .= "$chr\tRepeatMasker\tTransposon\t$start\t$end\t$score\t$strand\t.\tID=$id;Target=$target $target_start $target_end;Class=$class;PercDiv=$t[1];PercDel=$t[2];PercIns=$t[3];\n";
		$mark++;
	}
	close IN;

	open OUT,">$file.gff" || die "fail creat $file";
	print OUT "##gff-version 3\n$output";
	close OUT;	

}


##6.00e-22      117 WUBlastX Chr07frag1M       26352    26816    - RETRO1_pol      LTR/Gypsy            360      533
####################################################
sub annot_to_gff3 {
	my $file = shift;
	my $pre_tag = shift;
	my $output;

	$pre_tag .= "_" if($pre_tag); 
	my $line_num = `wc -l $file`;
	$line_num = $1 if($line_num =~ /^(\d+)/);
	my $mark = create_marker($line_num);

	open IN,$file || die "fail open $file";
	while (<IN>) {
		s/^\s+//;
		my @t = split /\s+/;
		next if($t[0] =~ /^pValue/);
		my $start = $t[4];
		my $end = $t[5];
		my $id = $pre_tag."TP".$mark;
		my $chr = $t[3];
		my $score = $t[1]; ##use pvalue here
		my $strand = $t[6];
		my $target = $t[7];
		my $class = $t[8];
		my @ary = ($t[9],$t[10]);
		@ary = sort {$a<=>$b} @ary;
		my ($target_start,$target_end) = ($ary[0],$ary[1]);

		$output .= "$chr\tRepeatProteinMask\tTEprotein\t$start\t$end\t$score\t$strand\t.\tID=$id;Target=$target $target_start $target_end;Class=$class;pValue=$t[0];\n";
		$mark++;
	}
	close IN;

	open OUT,">$file.gff" || die "fail creat $file";
	print OUT "##gff-version 3\n$output";
	close OUT;	

}
