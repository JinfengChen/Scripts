#!/usr/bin/perl

=head1 Name

stat_TE.pl  --  stat the TE content by type, subtype, and family

=head1 Description

This program can read repeatmasker .out file, RepeatProteinMask .annot file,
and gff file, to statistic the TE content for each type, subtype, and family.
Note that all the numbers are of non-redundant length.



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-12-21
  Note:

=head1 Usage
  
  perl stat_TE.pl [options]
  --repeat <file>   set repeatmasker result file
  --protein <file>  set RepeatProteinMask result file 
  --gff <file>      set the gff format input file
  --rank <str>      set the rank of stat: all, type, subtype and family. default= type
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

   perl ../bin/stat_TE.pl --repeat ./rice.frag1M.fa.RepeatMasker.out --protein rice.frag1M.fa.Proteinmask.annot

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my (@RepeatMasker,@RepeatProteinMask,@Gff,$Rank);
my ($Verbose,$Help);
GetOptions(
	"repeat:s"=>\@RepeatMasker,
	"protein:s"=>\@RepeatProteinMask,
	"gff:s"=>\@Gff,
	"rank:s"=>\$Rank,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Rank ||= "type";
die `pod2text $0` if ( $Help);

my %Data;

foreach my $file (@RepeatMasker) {
	Read_RepeatMasker($file,\%Data);
}

foreach my $file (@RepeatProteinMask) {
	Read_RepeatProteinMask($file,\%Data);
}

foreach my $file (@Gff) {
	Read_gff($file,\%Data);
}


foreach my $TE_type (sort keys %Data) {
	my $TE_type_p = $Data{$TE_type};
	my $TE_type_size = 0;
	foreach my $chr (sort keys %$TE_type_p) {
		$TE_type_size += Conjoin_fragment($TE_type_p->{$chr});
	}
	print "$TE_type\t$TE_type_size\n";
        #print "$TE_type\n";
}


####################################################
################### Sub Routines ###################
####################################################


##read RepeatMasker .out file
#usage: Read_RepeatMasker($file,\%hash);
############################################
sub Read_gff{
	my $file=shift;
	my $hash_p=shift; 
	
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		next if(/^\#/);
		my @temp=split(/\t/);

		my $tname = $temp[0];
		my $strand = $temp[6];
		my $start = $temp[3];
		my $end = $temp[4];
		#my $TE_name = $1 if ($temp[8] =~ /Target=([^;]+);*/);
                my $TE_name = $1 if ($temp[8] =~ /Target=(.*)\s+\d+\s+;*/);
		my $TE_class = $1 if ($temp[8] =~ /Class=([^;]+);*/);
		$TE_class =~ s/\?$//;
		
		

		my ($TE_type,$TE_subtype);
		$TE_type = $1 if($TE_class =~ /^([^\/]+)\/*/);
		$TE_subtype = $1 if($TE_class =~ /^[^\/]+\/([^\/]+)/);
		$TE_subtype ||= $TE_type;
		my $need_stat;
		if ($Rank eq "all") {
			$need_stat = "Total_TE";
		}elsif($Rank eq "type"){
			$need_stat = $TE_type;
		}elsif($Rank eq "subtype"){
			$need_stat = "$TE_type/$TE_subtype";
		}elsif($Rank eq "family"){
			$need_stat = "$TE_type/$TE_subtype/$TE_name";
		}

		push @{$hash_p->{$need_stat}{$tname}}, [$start,$end]; 
	}
	close(IN);

}


##read RepeatMasker .out file
#usage: Read_RepeatMasker($file,\%hash);
############################################
sub Read_RepeatMasker{
	my $file=shift;
	my $hash_p=shift; 
	
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		s/^\s+//;
		my @temp=split(/\s+/);
		
		next if($temp[8] ne 'C' && $temp[8] ne '+');
		my $tname = $temp[4];
		my $strand = ($temp[8] eq '+') ? '+' : '-';
		my $start = $temp[5];
		my $end = $temp[6];
		my $TE_name = $temp[9];
		my $TE_class = $temp[10];
		$TE_class =~ s/\?$//;
		
		my ($TE_type,$TE_subtype);
		$TE_type = $1 if($TE_class =~ /^([^\/]+)\/*/);
		$TE_subtype = $1 if($TE_class =~ /^[^\/]+\/([^\/]+)/);
		$TE_subtype ||= $TE_type;
		my $need_stat;
		if ($Rank eq "all") {
			$need_stat = "Total_TE";
		}elsif($Rank eq "type"){
			$need_stat = $TE_type;
		}elsif($Rank eq "subtype"){
			$need_stat = "$TE_type/$TE_subtype";
		}elsif($Rank eq "family"){
			$need_stat = "$TE_type/$TE_subtype/$TE_name";
		}

		push @{$hash_p->{$need_stat}{$tname}}, [$start,$end]; 
	}
	close(IN);

}



sub Read_RepeatProteinMask{
	my $file=shift;
	my $hash_p=shift; 
	
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		s/^\s+//;
		next if(/^pValue/);
		my @temp=split(/\s+/);
		my $tname = $temp[3];
		my $strand = $temp[6];
		next if($strand ne '-' && $strand ne '+');
		my $start = $temp[4];
		my $end = $temp[5];
		my $TE_name = $temp[7];
		my $TE_class = $temp[8];
		$TE_class =~ s/\?$//;
		
		my ($TE_type,$TE_subtype);
		$TE_type = $1 if($TE_class =~ /^([^\/]+)\/*/);
		$TE_subtype = $1 if($TE_class =~ /^[^\/]+\/([^\/]+)/);
		$TE_subtype ||= $TE_type;
		my $need_stat;
		if ($Rank eq "all") {
			$need_stat = "Total_TE";
		}elsif($Rank eq "type"){
			$need_stat = $TE_type;
		}elsif($Rank eq "subtype"){
			$need_stat = "$TE_type/$TE_subtype";
		}elsif($Rank eq "family"){
			$need_stat = "$TE_type/$TE_subtype/$TE_name";
		}

		push @{$hash_p->{$need_stat}{$tname}}, [$start,$end]; 
	}
	close(IN);

}


##conjoin the overlapped fragments, and caculate the redundant size
##usage: conjoin_fragment(\@pos);
##		 my ($all_size,$pure_size,$redunt_size) = conjoin_fragment(\@pos);
##Alert: changing the pointer's value can cause serious confusion.
sub Conjoin_fragment{
	my $pos_p = shift; ##point to the two dimension input array
	my $distance = shift || 0;
	my $new_p = [];         ##point to the two demension result array
	
	my ($all_size, $pure_size, $redunt_size) = (0,0,0); 
	
	return (0,0,0) unless(@$pos_p);

	foreach my $p (@$pos_p) {
			($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
			$all_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	@$pos_p = sort {$a->[0] <=>$b->[0]} @$pos_p;
	push @$new_p, (shift @$pos_p);
	
	foreach my $p (@$pos_p) {
			if ( ($p->[0] - $new_p->[-1][1]) <= $distance ) { # conjoin two neigbor fragements when their distance lower than 10bp
					if ($new_p->[-1][1] < $p->[1]) {
							$new_p->[-1][1] = $p->[1]; 
					}
					
			}else{  ## not conjoin
					push @$new_p, $p;
			}
	}
	@$pos_p = @$new_p;

	foreach my $p (@$pos_p) {
			$pure_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	$redunt_size = $all_size - $pure_size;
	return $pure_size;
}
