#!/usr/bin/perl

=head1 Name

ortholog_paralog_family.pl  -- construct gene families for mammalians

=head1 Description

This pipeline mainly use the Treefam's method and software, include the following steps:

1. an all-versus-all BLAST with proteins with Evalue cutoff (default 1e-7), note that each cds and protein pair must have the same gene_id in treefam's mode, such as BGIBMGA007261_BOMMO.

2. conjoin the blast alignments by solar, filter paired redundance in solar result, filter low alignments by size cutoff (default >= 1/3 to both genes), convert bit score to percent score (range 0 - 100) by algorithm: scoreP1P2 / max(scoreP1P1,scoreP2P2), and group the genes by hcluster_sg with default parameters(min_weight 10, min_density 0.34, max_size 500).

3. make gene family direcotry structure, using two-rank directory structure. start from 1, the format is like: 1-100, 101-200, 201-300, .......and put cds, pep, structure, and annotation files into each family.

4. do multiple sequence alignment using muscle, and covert pep alignment to cds alignment. We will use the -maxiters parameter to increase the speed for large family.

5. build gene phylogeny tree using treebest, which will execute nj-dn, nj-ds, nj-mm, phyml-aa, and phyml-nt algorithms, and merge these trees into a best consensus tree. It is also adjusted by the species tree.

6. calculate Ka/Ks, and 4dtv distance. Trace Ka/Ks onto each branch of the tree

7. draw figure include tree topology, multiple alignment, exon border, and annotation information.


=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 8.0,  Date: 2008-12-8

=head1 Usage
  
  $0 [options] <pep_file.fa> <cds_file.fa> <structure_file> <annotation_file> <category_file> <species_tree_file>
   --blast_evalue <str>   set the E-value of blast, default 1e-7
   --align_rate <num>     set the aligning rate for a pair of genes, default 0.33
   --min_weight <num>     set the minimum edge weight for hcluster_sg, default 10
   --min_density <num>    set the minimum edge density for hcluster_sg, default 0.34
   --max_size <num>       set the maximum family size for hcluster_sg, default 500
   --key_species <str>    set the key species, exa. BOMMO, default no
   --subdir_num <num>     set the number of subdir under gene_families dir, default 100
   --cpu <num>	          set the cpu number to use in parallel, default=3   
   --run <str>            set the parallel type, qsub, or multi, default=qsub
   --outdir <str>         set the result directory, default="./"
   --step <num>           set the start running step, default 1234567
   --synteny              generate synteny data, default no
   --distance             generate distance files inside and between species
   --single_copy          seperate the single copy nhx and figures, default no
   --verbose              output verbose information to screen  
   --help                 output help information to screen  

=head1 Exmple

  perl ../bin/ortholog_paralog_family.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species BOMMO --synteny --distance --single_copy &

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my ($blast_eval,$minimum_edge_weight,$minimum_edge_density,$maximum_size,$fam_in_subdir_num);
my ($Key_species,$Cpu,$Run,$Outdir,$Synteny,$Single_copy,$Distance,$Align_rate);
my ($Step);
my ($Verbose,$Help);
GetOptions(
	"blast_evalue:s"=>\$blast_eval,
	"align_rate:f"=>\$Align_rate,
	"min_weight:i"=>\$minimum_edge_weight,
	"min_density:f"=>\$minimum_edge_density,
	"max_size:i"=>\$maximum_size,
	"key_species:s"=>\$Key_species,
	"synteny"=>\$Synteny,
	"distance"=>\$Distance,
	"single_copy"=>\$Single_copy,
	"subdir_num:i"=>\$fam_in_subdir_num,
	"cpu:i"=>\$Cpu,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"step:i"=>\$Step,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Cpu ||= 4;
$Run ||= "multi";
$Outdir ||= "../output";
$Step ||= "1234567";
$blast_eval ||= 1e-7;  ##1e-10的错误率量级，可以保证1e7个比对结果都是正确的
$Align_rate ||= 0.33;
$minimum_edge_weight ||= 10; 
$minimum_edge_density ||= 0.34;
$maximum_size ||= 500;
$fam_in_subdir_num ||= 500; ##means 100 families in each subdirectory under top gene_families directory
die `pod2text $0` if (@ARGV < 6 || $Help);

my $pep_file = shift;
my $cds_file = shift;
my $structure_file = shift;
my $annotation_file = shift;
my $category_file = shift;
my $species_tree = shift; 

my $common_bin = "/home/jfchen/FFproject/tools/bin/";

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my $pep_file_basename = basename($pep_file);

my $formatdb = "/home/biosoftware/bin/formatdb";
my $blastall = "/home/biosoftware/bin/blastall";
my $solar = "/home/jfchen/FFproject/tools/solar/solar/solar.pl";
my $hcluster_sg = "/home/biosoftware/bin/hcluster_sg";
my $muscle = "/home/biosoftware/bin/muscle";
my $treebest = "/home/biosoftware/bin/treebest";
my $Kaks_Calculator = "/home/biosoftware/bin/kaks_calculator";


my $qsub_sge = "$common_bin/qsub-sge.pl";
my $multi_process = "$common_bin/multi-process.pl";
my $fastaDeal = "$common_bin/fastaDeal.pl";

my @fam; ##store all the families in two-demensional array
my %anno; ##store annotation information of each gene
my $fam_dir = "$Outdir/gene_families"; ##store the result families 

$Step = convert_step($Step) if($Step =~ /^\d-\d$/);

##step 1: running all vs all blast, using protein sequence, force to simplify the protein head lines
if ($Step =~ /1/) {
	`$formatdb -i $pep_file -p T -o T`;
	`perl $Bin/simplify_pep_fasta.pl $pep_file > $pep_file.simple; mv $pep_file.simple $pep_file`;
	`perl $fastaDeal -cutf $Cpu $pep_file -outdir $Outdir`;
	my @subfiles = glob("$Outdir/$pep_file_basename.cut/*.*");
	
	my $blast_shell = "$Outdir/blast_shell.sh";
	open OUT,">$blast_shell" || die "fail open $blast_shell\n";
	foreach my $subfile (@subfiles) {
		print OUT "$blastall -p blastp -m8 -e $blast_eval -F F -d $pep_file -i $subfile -o $subfile.blast.m8 \n";
	}
	close OUT;
	
	`perl $qsub_sge  --resource  --maxjob $Cpu -reqsub $blast_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $blast_shell` if ($Run eq "multi");
	
	`rm $pep_file.p??`;
	`rm $Outdir/all_vs_all.blast.m8` if(-e "$Outdir/all_vs_all.blast.m8");
	foreach my $subfile (@subfiles) {
		`cat $subfile.blast.m8 >> $Outdir/all_vs_all.blast.m8`;
	}

	print STDERR "all vs all protein blast done\n\n" if(defined $Verbose);
}


##step 2: cluster gene families
if ($Step =~ /2/) {
	my $category_genes_core = basename($category_file);
	my $category_genes_file = "$Outdir/$category_genes_core.genes";
	`perl $Bin/make_category_for_hcluster.pl $category_file $pep_file > $category_genes_file`;
	my $cluster_shell = "$Outdir/hcluster_shell.sh";
	open OUT,">$cluster_shell" || die "fail open $cluster_shell\n";

	print OUT "perl $solar -a prot2prot -f m8 -z $Outdir/all_vs_all.blast.m8 > $Outdir/all_vs_all.blast.m8.solar.raw; ",
		"perl $Bin/filter_solar_paired_redundance.pl $Outdir/all_vs_all.blast.m8.solar.raw > $Outdir/all_vs_all.blast.m8.solar.noPaired; ",
		"perl $Bin/filter_solar_align_rate.pl $pep_file $Outdir/all_vs_all.blast.m8.solar.noPaired $Align_rate > $Outdir/all_vs_all.blast.m8.solar; ",
		"perl $Bin/bitScore_to_hclusterScore.pl  $Outdir/all_vs_all.blast.m8.solar > $Outdir/all_vs_all.blast.m8.solar.forHC 2> $Outdir/all_vs_all.blast.m8.solar.forHC.warn; ",
		"perl $Bin/draw_score_distribution.pl $Outdir/all_vs_all.blast.m8.solar.forHC; ",
		"$hcluster_sg -w $minimum_edge_weight -s $minimum_edge_density -m $maximum_size -b 0.1 -C $category_genes_file $Outdir/all_vs_all.blast.m8.solar.forHC > $Outdir/all_vs_all.blast.m8.solar.forHC.hcluster ; \n";
	
	close OUT;

	`perl $qsub_sge  --maxjob $Cpu -reqsub $cluster_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $cluster_shell` if ($Run eq "multi");
	
	##filter and rename the cluster file(remove those family with single gene), and do some simple statistics
	my $hcluster_result_file = "$Outdir/all_vs_all.blast.m8.solar.forHC.hcluster";
	filter_hcluster($hcluster_result_file);
	`perl $Bin/hcluster_stat.pl $category_file $hcluster_result_file`;
	`perl $Bin/fishInWinter.pl  $hcluster_result_file.stat.single-copy $hcluster_result_file > $hcluster_result_file.single-copy`;
	
	##make pair wise list file for synteny analysis
	if (defined $Synteny) {
		
		my %species_list;
		my $synteny_dir = "$Outdir/synteny_data";
		my $synteny_shell = "$Outdir/synteny_shell.sh";
		read_category($category_file,\%species_list);
		open OUT,">$synteny_shell" || die "fail open $synteny_shell\n";
		foreach my $spec (sort keys %species_list) {
			print OUT "perl $Bin/make_pair_wise_list.pl $Outdir/all_vs_all.blast.m8.solar $structure_file $spec -outdir $synteny_dir\n";
		}
		foreach my $spec1 (sort keys %species_list) {
			foreach my $spec2 (sort keys %species_list) {
				next if($spec1 eq $spec2);
				print OUT "perl $Bin/make_pair_wise_list_twospecies.pl $Outdir/all_vs_all.blast.m8.solar $structure_file $spec1 $spec2 -outdir $synteny_dir\n";
			}
		}
		close OUT;
		
		`perl $qsub_sge  --maxjob $Cpu -reqsub $synteny_shell` if ($Run eq "qsub");
		`perl $multi_process -cpu $Cpu $synteny_shell` if ($Run eq "multi");
	}

	print STDERR "cluster gene families done\n\n" if(defined $Verbose);
}



##this function is needed by all the following steps

read_hcluster_family("$Outdir/all_vs_all.blast.m8.solar.forHC.hcluster",\@fam);


##step 3: make gene family direcotry structure, using two-rank directory structure
##start from 1, the format is like: 1-100, 101-200, 201-300, .......
if ($Step =~ /3/) {
	
	`mkdir $fam_dir` unless(-d $fam_dir);
	my $family_id = 1;
	foreach (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		mkpath("$fam_dir/$sub_dir/$family_id");
		$family_id++;
	}


	my %pep;
	Read_fasta($pep_file,\%pep);
	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		open OUT,">$fam_dir/$sub_dir/$family_id/$family_id.pep" || die "fail $family_id.pep";
		foreach my $gene_id (@$p) {
			print OUT ">",$pep{$gene_id}{head},"\n",$pep{$gene_id}{seq};	
		}
		close OUT;
		$family_id++;
	}
	undef %pep;

	my %cds;
	Read_fasta($cds_file,\%cds);
	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		open OUT,">$fam_dir/$sub_dir/$family_id/$family_id.cds" || die "fail $family_id.cds";
		foreach my $gene_id (@$p) {
			print OUT ">",$cds{$gene_id}{head},"\n",$cds{$gene_id}{seq};	
		}
		close OUT;
		$family_id++;
	}
	undef %cds;
	
	my %structure;
	Read_table($structure_file,\%structure);
	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		open OUT,">$fam_dir/$sub_dir/$family_id/$family_id.structure" || die "fail $family_id.structure";
		foreach my $gene_id (@$p) {
			print OUT "$gene_id\t$structure{$gene_id}\n";	
		}
		close OUT;
		$family_id++;
	}
	undef %structure;

	my %annotation;
	Read_table($annotation_file,\%annotation);

	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		open OUT,">$fam_dir/$sub_dir/$family_id/$family_id.annotation" || die "fail $family_id.annotation";
		foreach my $gene_id (@$p) {
			print OUT "$gene_id\t$annotation{$gene_id}\n";	
		}
		close OUT;
		$family_id++;
	}
	undef %annotation;
	print STDERR "gene family directory made\n\n" if(defined $Verbose);
}

##step 4: do multiple sequence alignment using muscle, and covert pep alignment to cds alignment
if ($Step =~ /4/) {
	my $family_id = 1;
	
	my $muscle_shell = "$Outdir/muscle_shell.sh";
	open OUT,">$muscle_shell" || die "fail open $muscle_shell\n";

	foreach my $p (@fam) {
		my $gene_num = scalar @$p;
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

		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		print OUT "$muscle $muscle_option -in $fam_dir/$sub_dir/$family_id/$family_id.pep -out $fam_dir/$sub_dir/$family_id/$family_id.pep.muscle; perl $Bin/pepMfa_to_cdsMfa.pl $fam_dir/$sub_dir/$family_id/$family_id.pep.muscle $fam_dir/$sub_dir/$family_id/$family_id.cds > $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle; \n";
		$family_id++;
	}
	
	close OUT;
	
	`perl $qsub_sge  --maxjob $Cpu -reqsub --lines 20 $muscle_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $muscle_shell` if ($Run eq "multi");

	print STDERR "muscle finished\n\n" if(defined $Verbose);
	
}

##step 5: build gene evolutionary tree using treebest with cds alignment
if ($Step =~ /5/) {
	my $family_id = 1;
	
	my $treebest_shell = "$Outdir/treebest_shell.sh";
	open OUT,">$treebest_shell" || die "fail open $treebest_shell\n";
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		print OUT "$treebest best -f $species_tree -p $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle -o $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx  $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle 2> $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.log; "; 
		print OUT "$treebest nj -c $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx -v -s $species_tree $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle > $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.ortho; ";
		print OUT "perl $Bin/get_ortholog_pairs.pl $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.ortho > $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.ortho.list; \n";
		$family_id++;
	}
	
	close OUT;
	
	`perl $qsub_sge  --maxjob $Cpu -reqsub --lines 20 $treebest_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $treebest_shell` if ($Run eq "multi");
	
	print STDERR "njtree finished\n\n" if(defined $Verbose);

}


##step 6: caculate Ka, ks, and ka/ks, 4dtv, and draw distribution figure for Ks and 4dtv
if ($Step =~ /6/) {
	my $kaks_shell = "$Outdir/kaks_to_tree.sh";
	open OUT,">$kaks_shell" || die "fail open $kaks_shell\n";

	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		print OUT "perl $Bin/add_kaks_to_tree.pl --treebest $treebest --kaks_caculator $Kaks_Calculator $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx;  ",
			"perl $Bin/calculate_4DTV_correction.pl $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt > $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt.4dtv;  ",
			"perl $Bin/calculate_cds_aa_identity.pl $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt > $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt.identity; \n";
		
		$family_id++;
	}
	close OUT;
	
	`perl $qsub_sge  --maxjob $Cpu -reqsub --lines 20 $kaks_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $kaks_shell` if ($Run eq "multi");
	
	print STDERR "kaks to tree finished\n\n" if(defined $Verbose);
	##generate the whole KaKs, 4dtv, and ortholog file, and draw distribution figures
	if (defined $Distance) {
		
		my $distance_dir = "$Outdir/distance_data";
		mkdir($distance_dir) unless(-d $distance_dir);

		my $family_id = 1;
		foreach my $p (@fam) {
			my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
			`cat $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt.4dtv | awk '\$4>0 && !/^tag/' >> $distance_dir/all.4dtv`;
			`cat $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt.KaKs | awk '\$4!="Ks" && \$4!="nan" && \$4!="NA"' >> $distance_dir/all.KaKs`;		
			`cat $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.axt.identity | awk '!/^tag/' >> $distance_dir/all.identity`;
			`cat $fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.ortho.list >> $distance_dir/all.ortho`;
			$family_id++;
		}
	}

}



##step 7: draw svg figure for tree, phylogeny tree, alignment structure (exon-intron border), text description
if ($Step =~ /7/) {
	
	my $draw_svg_shell = "$Outdir/draw_svg_tree.sh";
	open OUT,">$draw_svg_shell" || die "fail open $draw_svg_shell\n";

	my $family_id = 1;
	foreach my $p (@fam) {
		my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
		my $tree_nhx_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.nhx";
		my $annotation_file = "$fam_dir/$sub_dir/$family_id/$family_id.annotation";
		my $structure_file = "$fam_dir/$sub_dir/$family_id/$family_id.structure";
		my $cds_muscle_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle";
		my $tree_id_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.list";
		my $tree_svg_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.svg";

		##we will develop another program to show the alignment and annotation information
		if(-f $tree_nhx_file){
			print OUT "perl $Bin/draw_tree_svg.pl $tree_nhx_file $Key_species; ",
				"perl  $Bin/display_muscle_ouput.pl  --seq_file $cds_muscle_file --annotation_file $annotation_file  --structure_file $structure_file --id_file $tree_id_file --exon_mark line --figure_width 1000 --left_margin 510 --right_margin 30  --top_margin 30 --bottom_margin 30 --join T   --coordinate --outdir $fam_dir/$sub_dir/$family_id; ",
				"perl $Bin/combine_two_svg.pl $tree_svg_file $cds_muscle_file.svg > $tree_svg_file.info.svg;\n";
		}else{
			warn "$tree_nhx_file not exist!"
		}
		$family_id++;
	}
	close OUT;

	`perl $qsub_sge  --maxjob $Cpu -reqsub --lines 20 $draw_svg_shell` if ($Run eq "qsub");
	`perl $multi_process -cpu $Cpu $draw_svg_shell` if ($Run eq "multi");
	
	##separate the single-copy figures
	if (defined $Single_copy) {
		my %single_copy_list;
		read_single_copy("$Outdir/all_vs_all.blast.m8.solar.forHC.hcluster.stat.single-copy",\%single_copy_list);
		my $single_figure_dir = "$Outdir/single_copy";
		mkdir($single_figure_dir) unless(-d $single_figure_dir);
		my $family_id = 1;
		foreach my $p (@fam) {
			my $sub_dir = cacul_sub_dir($family_id,$fam_in_subdir_num);
			my $svg_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.svg";
			my $nhx_file = "$fam_dir/$sub_dir/$family_id/$family_id.cds.muscle.best.nhx.Ks.Ka.KaKs.nhx";
			`cp $svg_file $nhx_file $single_figure_dir` if(exists $single_copy_list{$family_id});
			$family_id++;
		}
	}

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
		$hash->{$t[0]} = $t[1] if($t[0] && $t[1]);
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
		$hash->{$t[0]} = $t[1] if($t[0] && $t[1]);
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
