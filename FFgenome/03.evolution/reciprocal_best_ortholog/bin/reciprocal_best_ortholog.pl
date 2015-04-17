#!/usr/bin/perl

=head1 Name

reciprocal_best_ortholog.pl  -- use the reciprocal best method to identify ortholog genes

=head1 Description

This program is used to compare two species using the proteins, by finding the reciprocal
best hits. It will statistic how much of the proteins are aligned to each species, and caculate
the mean(median) identity. 

Besides the reciprocal best hit file, this program will also generate a file containing all the 
hits. Both in the standard format, which can be used as input for synteny analysis programs.

The input are two protein files in fasta format, each represent the whole proteome for a species.
The head lines must be like the following, with attribute [locus=chromosome:start:end:strand]
>Scaffold45_orf00004  locus=Scaffold45:2123:2965:-


=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 2.0,  Date: 2008-12-19

=head1 Usage
  
  perl reciprocal_best_ortholog.pl <pep1.fa>  <pep2.fa> 
  --eval <str>    set e-value cutoff for blast hits, default=1e-5
  --step <str>    start point of program running, default=1
  --outdir <str>  set result directory, default=./
  --verbose       output verbose information to screen  
  --help          output help information to screen  

=head1 Exmple

perl ../bin/reciprocal_best_ortholog.pl ../input/NC_003155.gb.pep   ../input/NC_003888.gb.pep 

=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
#use lib "$Bin/../../../common_bin";
#use GACP qw(parse_config);

my ($Cutoff,$Step,$Outdir,$Verbose,$Help);
GetOptions(
	"eval:s"=>\$Cutoff,
	"step:i"=>\$Step,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Step ||= 1;
$Cutoff ||= 1e-5;
$Outdir =~ s/\/$//;
$Outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $pep1_file = shift;
my $pep2_file = shift;

my $pep1_file_basename = basename($pep1_file);
my $pep2_file_basename = basename($pep2_file);

my (%Best_pair,%best_pair1,%best_pair2,$Best_pair_num);
my (%Seq1_pos,$Seq1_gene_num,%Seq2_pos,$Seq2_gene_num);

##set the path of softwares
#my $config_file = "$Bin/../../../config.txt";
#my $common_bin = "$Bin/../../../common_bin";
my $common_bin = $Bin;
#my $blastall = parse_config($config_file,"blastall")." -b 10000 -v 10000 -F F";;
#my $formatdb = parse_config($config_file,"formatdb");
my $blastall = "blastall -F F";
my $formatdb = "formatdb";
##run blastp between two protein sets, this step consume most of the time
if ($Step == 1) {
	`$formatdb -i $pep1_file -p T & $formatdb -i $pep2_file -p T & wait`; 
	`$blastall -p blastp -a 20 -e $Cutoff -m 8 -d $pep2_file -i $pep1_file -o $Outdir/$pep1_file_basename--$pep2_file_basename.m8 & $blastall -p blastp -a 20 -e $Cutoff -m 8 -d $pep1_file -i $pep2_file -o $Outdir/$pep2_file_basename--$pep1_file_basename.m8 & wait`;
	#`rm $pep1_file.p??  $pep2_file.p??`;
}

##find reciprocal best hit pair and do some statistics
if ($Step <= 2) {

	##get the best hit with the lowest E-value
	`perl $Bin/bestAlign.pl -f m8 -cutoff $Cutoff $Outdir/$pep1_file_basename--$pep2_file_basename.m8 > $Outdir/$pep1_file_basename--$pep2_file_basename.m8.best`;
	`perl $Bin/bestAlign.pl -f m8 -cutoff $Cutoff $Outdir/$pep2_file_basename--$pep1_file_basename.m8 > $Outdir/$pep2_file_basename--$pep1_file_basename.m8.best`;
	
	##get the reciprocal best hit ortholog gene pairs
	read_m8_best_pair("$Outdir/$pep1_file_basename--$pep2_file_basename.m8.best",\%best_pair1);
	read_m8_best_pair("$Outdir/$pep2_file_basename--$pep1_file_basename.m8.best",\%best_pair2);
	foreach my $pep1 (sort keys %best_pair1) {
		my $pep2 = $best_pair1{$pep1};
		if (exists $best_pair2{$pep2} && $best_pair2{$pep2} eq $pep1) {
			$Best_pair{$pep1} = $pep2;
			$Best_pair_num++;
		}
	}
	#print Dumper \%Best_pair;
	
	Read_pep_fasta($pep1_file,\%Seq1_pos,\$Seq1_gene_num);
	Read_pep_fasta($pep2_file,\%Seq2_pos,\$Seq2_gene_num);
	
	my (%identity,%evalue);
	my ($mean_identity,$median_identity)= read_m8_best_identity("$Outdir/$pep1_file_basename--$pep2_file_basename.m8.best",\%identity);
	read_m8_best_eval("$Outdir/$pep1_file_basename--$pep2_file_basename.m8.best",\%evalue);

	my $Seq1_gene_num_rate = sprintf("%.2f",$Best_pair_num / $Seq1_gene_num * 100);
	my $Seq2_gene_num_rate = sprintf("%.2f",$Best_pair_num / $Seq2_gene_num * 100);
	
	##generate the statistic file for ortholog genes
	my $reciprocal_best_stat = "$Outdir/".$pep1_file_basename.".".$pep2_file_basename.".reciprocal.best.stat";
	open OUT, ">$reciprocal_best_stat" || die "fail creat $reciprocal_best_stat";
	print OUT 
		"Reciprocal best hits, Evalue_cutoff: $Cutoff\n",
		"Query:  $Seq1_gene_num genes, $Best_pair_num ($Seq1_gene_num_rate\%) aligned\n",
		"Target: $Seq2_gene_num genes, $Best_pair_num ($Seq2_gene_num_rate\%) aligned\n",
		"identity(mean median):  $mean_identity  $median_identity\n";
	
	##generate the position list file for ortholog genes  
	my $reciprocal_best_list = "$Outdir/".$pep1_file_basename.".".$pep2_file_basename.".reciprocal.best.list";
	open OUT, ">$reciprocal_best_list" || die "fail creat $reciprocal_best_list";
	foreach my $gene1_id (sort keys %Best_pair) {
		my $gene2_id = $Best_pair{$gene1_id};

		my $ide = $identity{$gene1_id}{$gene2_id};
		my $eva = $evalue{$gene1_id}{$gene2_id};
		print OUT 
			"$gene1_id\t$Seq1_pos{$gene1_id}[3]\t$Seq1_pos{$gene1_id}[2]\t$Seq1_pos{$gene1_id}[0]\t$Seq1_pos{$gene1_id}[1]\t",
			"$gene2_id\t$Seq2_pos{$gene2_id}[3]\t$Seq2_pos{$gene2_id}[2]\t$Seq2_pos{$gene2_id}[0]\t$Seq2_pos{$gene2_id}[1]\t",
			"$eva\t$ide\n";
	}
	
	close OUT;


	##generate the all hit file for synteny analysis
	##order: [$gene_start,$gene_end,$strand,$chr]
	my $all_hits_list = "$Outdir/".$pep1_file_basename.".".$pep2_file_basename.".all.hits.list";
	open OUT, ">$all_hits_list" || die "fail open $all_hits_list";
	open IN,"$Outdir/$pep1_file_basename--$pep2_file_basename.m8" || die "fail open";
	while (<IN>) {
		my @t = split /\t/;
		my $gene1_id = $t[0];
		my $gene2_id = $t[1];
		next if($gene1_id eq $gene2_id);
		print OUT 
			"$gene1_id\t$Seq1_pos{$gene1_id}[3]\t$Seq1_pos{$gene1_id}[2]\t$Seq1_pos{$gene1_id}[0]\t$Seq1_pos{$gene1_id}[1]\t",
			"$gene2_id\t$Seq2_pos{$gene2_id}[3]\t$Seq2_pos{$gene2_id}[2]\t$Seq2_pos{$gene2_id}[0]\t$Seq2_pos{$gene2_id}[1]\t",
			"$t[10]\t$t[2]\n";
                 
	}
	close IN;
	close OUT;
	
	##generate protein identity distribution figure
	 #`awk '{print \$12}' $reciprocal_best_list > $reciprocal_best_list.identity`;
	 #`perl $common_bin/distribute_pre.pl --frequence  --binsize 1  --header rect --color blue --Note "Distribution of protein identity"  --Y "Number of orthologs" --X "Protein identity (window=1)" $reciprocal_best_list.identity > $reciprocal_best_list.identity.lst`;
	 #`perl $common_bin/distribute_svg.pl $reciprocal_best_list.identity.lst $reciprocal_best_list.identity.svg`;
}
	



####################################################
################### Sub Routines ###################
####################################################


sub mean_median{
	my $ary_p = shift;
	my ($number,$total,$mean,$median);
	foreach  (@$ary_p) {
		$total += $_;
		$number++;
	}
	$mean=sprintf("%.2f",$total/$number);
	$median=sprintf("%.2f",$ary_p->[int $number/2] );
	
	return ($mean,$median);
}


##The order of fields for BLAST result in tabular format is: query id, database sequence (subject) id, percent identity, alignment length, number of mismatches, number of gap openings, query start, query end, subject start, subject end, Expect value, HSP bit score.
sub read_m8_best_identity {
	my $file = shift;
	my $hash_p = shift;

	my @val;
	open IN,$file || die "fail $file";
	while (<IN>) {
		my @t = split /\t/;
		$hash_p->{$t[0]}{$t[1]} = $t[2];
		push @val,$t[2];
	}
	close IN;

	my ($mean,$median) = mean_median(\@val);
	$mean = sprintf("%.2f",$mean);
	$median = sprintf("%.2f",$median);
	return ($mean,$median);
}


sub read_m8_best_eval {
	my $file = shift;
	my $hash_p = shift;

	my @val;
	open IN,$file || die "fail $file";
	while (<IN>) {
		my @t = split /\t/;
		$hash_p->{$t[0]}{$t[1]} = $t[10];
		push @val,$t[10];
	}
	close IN;

	my ($mean,$median) = mean_median(\@val);
	$mean = sprintf("%.2e",$mean);
	$median = sprintf("%.2e",$median);
	return "$mean  $median";
}


sub read_m8_best_pair {
	my $file = shift;
	my $hash_p = shift;
	open IN,$file || die "fail $file";
	while (<IN>) {
		my @t = split /\t/;
		$hash_p->{$t[0]} = $t[1];
	}
	close IN;
}


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_pep_fasta{
	my $file=shift;
	my $hash_p=shift;
	my $total_p = shift;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $head = $_;
		my ($gene_id,$chr,$gene_start,$gene_end,$strand);
		$gene_id = $1 if($head =~ /^(\S+)/);
		if($head =~ /([^:=]+):(\d+):(\d+):([+-])/){
			($chr,$gene_start,$gene_end,$strand) = ($1,$2,$3,$4);
		}else{
			($chr,$gene_start,$gene_end,$strand) = ("NO","NO","NO","NO");
		}

		$hash_p->{$gene_id} = [$gene_start,$gene_end,$strand,$chr];
		
		$/=">";
		my $seq = <IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";

		$$total_p++;
	}
	close(IN);
	
}
