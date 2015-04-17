#!/usr/bin/perl

=head1 Name

getGene.pl  --  get gene elements from sequences according to coordinates

=head1 Description

This program is used to get out the gene elements, such as exon, intron, mrna, gene, splicing site, 
5'-flanking, 3'-flanking. Both cds and gene can be used in combination with 5'-flanking and 3'-flanking.

The position file can be psl or gff format now. The genome sequence file must
be fasta format. Note that this program is originally designed for CDS region, in other words,
not consider the UTR regions. You should alter on this, in order not to make mistake. 
For psl format, it only considers about blocks. For gff format, it only recognize mRNA and CDS feature.
The program will detect the file format automatically, or you can set by "--posformat" manually.

The default for "--type" option is "mrna", here it has the same meaning as cds, because we don't consider
UTR region.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2007-12-26

=head1 Usage
  % $0  [option] <pos_file> <seq_file>
  --posformat <str>   specify position file format
  --type <str>        specify element type: exon,intron,mrna,gene,splice,flank5,flank3, default=mrna
  --flank3 <num>      specify 3'-flanking length, co-used with gene/mrna
  --flank5 <num>      specify 5'-flanking length, co-used with gene/mrna
  --verbose           output verbose information to screen  
  --help              output help information to screen  

=head1 Exmple

 perl ./getGene.pl chr01.psl chr01.fa 
 perl ./getGene.pl chr01.gff chr01.fa -type splice
 perl ./getGene.pl chr01.psl chr01.fa -type intron
 perl ./getGene.pl chr01.psl chr01.fa -type exon
 perl ./getGene.pl chr01.psl chr01.fa -type gene
 perl ./getGene.pl chr01.psl chr01.fa -type mrna -flank5 50 -flank3 50

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Posformat,$Type,$Flank5,$Flank3);
my ($Verbose,$Help);
GetOptions(
	"type:s"=>\$Type,
	"flank5:i"=>\$Flank5,
	"flank3:i"=>\$Flank3,
	"posformat:s"=>\$Posformat,
	"verbose!"=>\$Verbose,
	"help!"=>\$Help
);
$Type ||= "mrna";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $pos_file = shift;
my $seq_file = shift;

my %gene;

read_psl($pos_file,\%gene) if($Posformat eq 'psl' || $pos_file =~ /.psl$/);
read_gff($pos_file,\%gene) if($Posformat eq 'gff' || $pos_file =~ /.gff$/);

open(IN,$seq_file)||die("failed $seq_file\n");

$/=">"; <IN>; $/="\n";	
while (<IN>) {
	my $output;
	my $chr=$1 if(/^(\S+)/);
        chomp $chr;
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";
	my $seq_len=length($seq);
        $chr =$1 if ($chr=~/\_(chr\d+)\_/);
        print "$chr\t$seq_len\n";
	my $chr_pp=$gene{$chr};
	foreach  my $gene (sort keys %$chr_pp) {
                print "$gene\n";
		my $strand=$$chr_pp{$gene}{strand};
		next if(!exists $chr_pp->{$gene}{exon});
		my @exon = @{$chr_pp->{$gene}{exon}};
		
		if ($Type eq 'flank5') {
			my ($left_leng, $right_leng,$flank_str);
			if ($strand eq '+') {
				$left_leng = $Flank5 if ($strand eq '+' && $Flank5 );
				$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);
				$flank_str = substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
			}
			if ($strand eq '-') {
				$right_leng = $Flank5 if ($strand eq '-' && $Flank5 );
				$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);
				$flank_str = substr($seq,$exon[-1][1],$right_leng) if($right_leng);
				$flank_str = Complement_Reverse($flank_str);
			}
			Display_seq(\$flank_str);
			$output .= ">".$gene."   [flank5:$Flank5]\n".$flank_str;
		}

		if ($Type eq 'flank3') {
			my ($left_leng, $right_leng,$flank_str);
			if ($strand eq '-') {
				$left_leng = $Flank3 if ($strand eq '-' && $Flank3 );
				$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);
				$flank_str = substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
				$flank_str = Complement_Reverse($flank_str);
			}
			if ($strand eq '+') {
				$right_leng = $Flank3 if ($strand eq '+' && $Flank3 );
				$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);
				$flank_str = substr($seq,$exon[-1][1],$right_leng) if($right_leng);
				
			}
			Display_seq(\$flank_str);
			$output .= ">".$gene."   [flank3:$Flank3]\n".$flank_str;
		}

		
		if ($Type eq "exon") {
			my $mark;
			for (my $i=0; $i<@exon; $i++) {
				my $exon = substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);
				$exon = Complement_Reverse($exon) if($strand eq '-');
				$mark = $gene."_E".($i+1) if($strand eq '+');
				$mark = $gene."_E".(@exon - $i) if($strand eq '-');
				Display_seq(\$exon);
				$output .= ">".$mark."\n".$exon;
			}
		}

		if ($Type eq "mrna") {
			my $mrna;
			my ($left_leng, $right_leng);
			$left_leng = $Flank5 if ($strand eq '+' && $Flank5 );
			$left_leng = $Flank3 if ($strand eq '-' && $Flank3 );
			$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);

			$right_leng = $Flank3 if ($strand eq '+' && $Flank3 );
			$right_leng = $Flank5 if ($strand eq '-' && $Flank5 );
			$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);

			$mrna .= substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
			for (my $i=0; $i<@exon; $i++) {
				$mrna .= substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);	
			}
			$mrna .= substr($seq,$exon[-1][1],$right_leng) if($right_leng);

			$mrna = Complement_Reverse($mrna) if($strand eq '-');
			Display_seq(\$mrna);
			my $mark = "$gene  [mRNA]";
			$output .= ">".$mark."\n".$mrna;
		}

		if ($Type eq "gene") {
			my $geneseq;
			my ($left_leng, $right_leng);
			$left_leng = $Flank5 if ($strand eq '+' && $Flank5 );
			$left_leng = $Flank3 if ($strand eq '-' && $Flank3 );
			$left_leng = $exon[0][0]-1 if($left_leng > $exon[0][0]-1);

			$right_leng = $Flank3 if ($strand eq '+' && $Flank3 );
			$right_leng = $Flank5 if ($strand eq '-' && $Flank5 );
			$right_leng = $seq_len - $exon[-1][1] if($right_leng > $seq_len - $exon[-1][1]);

			$geneseq .= substr($seq,$exon[0][0]-$left_leng-1,$left_leng) if($left_leng);
			$geneseq .= substr($seq,$exon[0][0] - 1, $exon[-1][1] - $exon[0][0] + 1);
			$geneseq .= substr($seq,$exon[-1][1],$right_leng) if($right_leng);
			
			$geneseq = Complement_Reverse($geneseq) if($strand eq '-');
			Display_seq(\$geneseq);
			my $mark = "$gene  [gene]";
			$output .= ">".$mark."\n".$geneseq;
		}
		
		if ($Type eq "intron") {
			my $mark;
			for (my $i=0; $i<@exon-1; $i++) {
				my $intron = substr($seq,$exon[$i][1], $exon[$i+1][0]-$exon[$i][1]-1);
				$intron = Complement_Reverse($intron) if($strand eq '-');
				$mark = $gene."_I".($i+1) if($strand eq '+');
				$mark = $gene."_I".(@exon-1-$i) if($strand eq '-');
				Display_seq(\$intron);
				$output .= ">".$mark."\n".$intron;
			}
		}

		if ($Type eq "splice") {
			my $mark;
			for (my $i=0; $i<@exon-1; $i++) {
				my $intron = substr($seq,$exon[$i][1], $exon[$i+1][0]-$exon[$i][1]-1);
				$intron = Complement_Reverse($intron) if($strand eq '-');
				$mark = $gene."_I".($i+1) if($strand eq '+');
				$mark = $gene."_I".(scalar(@exon)-1-$i) if($strand eq '-');	
				my $splice = substr($intron,0,2)."/".substr($intron,length($intron)-2,2);
				$output .= $gene."\t".$mark."\t".$splice."\n";
			}
		}

		
		
		
	}
	print $output;
}

close(IN);




####################################################
################### Sub Routines ###################
####################################################

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################


#############################################
sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;

}
#############################################



sub read_psl{
	my $file=shift;
	my $ref=shift;
	open(REF,$file)||die("fail to open $file\n");
	while (<REF>) {
		my @temp=split(/\s+/,$_);
		my $chr=$temp[13];
		my $gene=$temp[9];
		my $strand=$temp[8];
		my @sizes=split(/,/,$temp[18]);
		my @starts=split(/,/,$temp[20]);
		my @exon;
		for (my $i=0; $i<@starts; $i++) {
			push @exon, [ $starts[$i]+1, $starts[$i]+$sizes[$i] ];
		}

		$$ref{$chr}{$gene}{strand}=$strand;
		$$ref{$chr}{$gene}{exon}=\@exon;


	}
	close(REF);
}

sub read_gff{
	my $file=shift;
	my $ref=shift;
	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		s/^\s+//;
		s/\s+$//;
		my @t = split(/\t/);
		my $tname = $t[0];
		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS') {
			$qname = $1 if($t[8] =~ /^GenePrediction\s+(\S+)/ || $t[8] =~ /^ID=([^;]+);*/ || $t[8] =~ /^Parent=([^;]+);*/);
			#print $qname."\n";
                        #$qname = $1 if($t[8] =~ /^ID=(.*?);/ || $t[8] =~ /Parent=(.*?);/);
		}
		if ($t[2] eq 'match' || $t[2] eq 'HSP') {
			$qname = $1 if($t[8] =~ /Target\s+\"(\S+)\"/);
		}

		
		if ($t[2] eq 'mRNA' || $t[2] eq 'match') {
			$ref->{$tname}{$qname}{strand} = $t[6];
		}
		if ($t[2] eq 'CDS' || $t[2] eq 'HSP') {
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]];
		}
	}
	close(IN);

	##print Dumper $ref;
	
	##change the exon order
	foreach my $chr (keys %$ref) {
		my $chr_p = $ref->{$chr};
		foreach my $gene (keys %$chr_p) {
                        #print "$gene\n";
			my $gene_p = $chr_p->{$gene};
			next if(!exists $gene_p->{exon});
			my @exon = sort {$a->[0] <=> $b->[0]} @{$gene_p->{exon}};
			$gene_p->{exon} = \@exon;
		}
	}
}
