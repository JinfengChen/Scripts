#!/usr/bin/perl

=head1 Name

ortholog_paralog_family_single.pl  -- build phylogeny tree and infer orthologs and paralogs for a single family

=head1 Description

This pipeline should be co-used with ortholog_paralog_family.pl; It works when you had changed the data 
of some families, and need to re-run the later analysis.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-10-10

=head1 Usage
  
  $0 [options]  <cds_file.fa> <structure_file> <annotation_file>  <species_tree_file>
   --key_species <str>    set the key species, exa. BOMMO, default no
   --outdir <str>    set the result directory, default="./"
   --step <num>      set the start running step, default 1234
   --verbose         output verbose information to screen  
   --help            output help information to screen  

=head1 Exmple

  perl ../../../../bin/ortholog_paralog_family_single.pl 1.cds 1.structure 1.annotation ../../../../input/species.nh --key_species BOMMO &

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory
use lib "$Bin/../../../common_bin";
use GACP qw(parse_config);

my ($Key_species,$Outdir,$Step);
my ($Verbose,$Help);
GetOptions(
	"key_species:s"=>\$Key_species,
	"outdir:s"=>\$Outdir,
	"step:i"=>\$Step,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Step ||= "1234";
die `pod2text $0` if (@ARGV < 4 || $Help);

my $cds_file = shift;
my $structure_file = shift;
my $annotation_file = shift;
my $species_tree = shift; 
my $cds_file_core = basename($cds_file);

$Outdir ||= dirname($cds_file);
$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my $pep_file = "$Outdir/$cds_file_core";
$pep_file =~ s/\.cds$/\.pep/;
my $cds_muscle_file = "$Outdir/$cds_file_core.muscle";
my $best_nhx_file = "$cds_muscle_file.best.nhx";

my $config_file = "$Bin/../../../config.txt";

my $muscle = parse_config($config_file,"muscle");
my $treebest = parse_config($config_file,"treebest");
my $Kaks_Calculator = parse_config($config_file,"kaks_caculator");

$Step = convert_step($Step) if($Step =~ /^\d-\d$/);


##step 1: do multiple sequence alignment using muscle
## and covert pep alignment to cds alignment
if ($Step =~ /1/) {
	
	`perl $Bin/cds2aa.pl $cds_file > $pep_file`;
	
	my $gene_num = `grep -c '>' $pep_file`;
	chomp $gene_num;
	my $muscle_option;
	if ($gene_num <= 20) {
		$muscle_option = '';
	}elsif ($gene_num > 20 && $gene_num < 50) {
		$muscle_option = "-maxiters 8";
	} elsif ($gene_num >= 50 && $gene_num < 200) {
		$muscle_option = "-maxiters 4";
	} elsif ($gene_num >= 200) {
		$muscle_option = "-maxiters 2";
	}
	`$muscle $muscle_option -in $pep_file -out $pep_file.muscle`; 
	`perl $Bin/pepMfa_to_cdsMfa.pl $pep_file.muscle $cds_file > $cds_muscle_file`;
	
	print STDERR "muscle finished\n\n" if(defined $Verbose);
	
}



##step 2: build gene evolutionary tree using njtree best with cds alignment
if ($Step =~ /2/) {
	
		`$treebest best -f $species_tree -p $cds_muscle_file -o $best_nhx_file $cds_muscle_file 2> $cds_muscle_file.best.log`; 
		
		`$treebest nj -c $best_nhx_file  -v -s $species_tree $cds_muscle_file > $best_nhx_file.ortho`;
		
		`perl $Bin/get_ortholog_pairs.pl $best_nhx_file.ortho > $best_nhx_file.ortho.list`;
		
	print STDERR "njtree finished\n\n" if(defined $Verbose);

}



##step 3: caculate Ka, ks, and ka/ks, 4dtv, and draw distribution figure for Ks and 4dtv
if ($Step =~ /3/) {
	

	`perl $Bin/add_kaks_to_tree.pl --treebest $treebest --kaks_caculator $Kaks_Calculator $cds_muscle_file $best_nhx_file`;

	`perl $Bin/calculate_4DTV_correction.pl $cds_muscle_file.axt > $cds_muscle_file.axt.4dtv`;

	print STDERR "kaks to tree finished\n\n" if(defined $Verbose);

}



##step 4: draw svg figure for tree
if ($Step =~ /4/) {
	
	`perl $Bin/draw_tree_svg.pl $best_nhx_file.Ks.Ka.KaKs.nhx $annotation_file $structure_file $Key_species`;

	print STDERR "tree drawing finished\n" if(defined $Verbose);
}




####################################################
################### Sub Routines ###################
####################################################

sub read_single_copy {
	my $file = shift;
	my $hash = shift;

	open IN,$file || die "fail $file";
	while (<IN>) {
		chomp;
		my @t = split /\s+/;
		$hash->{$t[0]} = 1;
	}
	close IN;
}

##read the species from category file
sub read_category{
	my $file = shift;
	my $hash = shift;

	open IN,$file || die "fail";
	while (<IN>) {
		chomp;
		my @t = split /\s+/;
		$hash->{$t[0]} = $t[1];
	}
	close IN;
}


##convert format "2-5" to format "2345"
#############################################
sub convert_step {
	my $str = shift;
	my $result;
	my ($start,$end) = ($1,$2) if($str =~ /^(\d)-(\d)$/);
	for (my $i=$start; $i<=$end; $i++) {
		$result .= "$i";
	}
	return $result;
}


#############################################
sub Read_table{
	my $file = shift;
	my $hash = shift;

	open IN,$file || die "fail open $file";
	while (<IN>) {
		chomp;
		my @t = split /\t/;
		$hash->{$t[0]} = $t[1];
	}
	close IN;
}


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



##caculate subdir id in a two ranks direcotry
##the format is: 1-100, 101-200, 201,300, .......
####################################################
sub cacul_sub_dir {
	my $family_id = shift;
	my $fam_in_subdir_num = shift;
	
	$family_id -= 1;
	my $sub_dir_num = ($family_id - $family_id % $fam_in_subdir_num) / $fam_in_subdir_num;
	my $start_fam_id = $sub_dir_num * $fam_in_subdir_num + 1;
	my $end_fam_id = ($sub_dir_num + 1) * $fam_in_subdir_num;

	my $sub_dir = "$start_fam_id-$end_fam_id";
	return $sub_dir;
}

sub read_hcluster_family{
	my $file = shift;
	my $ary_p = shift;

	open IN,$file || die "fail open $file";
	while (<IN>) {
		chomp;
		my @t = split /\t/;
		$t[-1] =~ s/,$//;
		my @k = split /,/, $t[-1];
		push @$ary_p,\@k;
	}
	close IN;
	
}

sub filter_hcluster{
	my $file = shift;
	
	my $family_id = 1;
	##filter and rename the cluster ID
	open OUT,">$file.filter" || die "fail $file.filter";
 	open IN,$file || die "fail open $file";
	while (<IN>) {
		chomp;
		my @t = split /\t/;
		next if($t[5] <= 1);
		$t[0] = $family_id;
		print OUT join("\t",@t)."\n";
		$family_id++;
	}
	close IN;
	close OUT;

	`mv $file.filter $file`;
}

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
		#$seq=~s/\s//g;
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


#read fasta file
#usage: Read_fasta($file,\%hash);
#############################################
sub Read_fasta_anno{
	my $file=shift;
	my $hash_p=shift;
	
	my $total_num;
	open(IN, $file) || die ("can not open $file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my ($name,$anno);
		$name = $1 if(/^(\S+)/);
		$anno = $2 if(/^(\S+)\s+(.+)/);
		$hash_p->{$name} = $anno;

		$/=">";
		my $seq = <IN>;
		#chomp $seq;
		#$seq=~s/\s//g;
		$/="\n";

		$total_num++;
	}
	close(IN);
	
	return $total_num;
}
