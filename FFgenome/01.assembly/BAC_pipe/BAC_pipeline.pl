#!user/bin/perl

=head1 Description
	To use the BAC sequence to valuate the quality of de novo assembly. Based on the pipeline of Cai Qingle, add the indel program from lujianliang.

=head1 Contact and Version
	Contact: liubinghang@genomics.org.cn
	Version: 2.0
	Date: 2010.2.21

=head1 Usage
  perl BAC_pepline.pl [options]
  --len1 <file>         set the scaffold length file
  --len2 <file>         set the bac length file
  --seq1 <file>         set the scaffold sequence file
  --seq2 <file>         set the bac sequence file
  --step <int>	1:prepare; 2:compare bac.seq with scaf.seq quickly,3:get pair list,caculate coverage,4:get pair seq ,blast, solar,get statistics;5. draw pictures  
  --gap1 <file>         draw contigs for one scaffold set (for step 5)
  --gap2 <file>         draw contigs for bac  set (for step 5)
  --repeat	<file>  set the repeat information of the bac sequence (for step 5)
  --gene	<file>  set the gene information of the bac sequence (for step 5)
  --outdir <str>        set the output directory
  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

perl /panfs/GAG/liubinghang/package/qc_assembly/BAC/BAC_pipeline.pl --seq1 /panfs/GAG/liubinghang/qc_assembly/crab_eating/5
00bac/test/bac.fa --seq2 /panfs/GAG/liubinghang/qc_assembly/crab_eating/500bac/test/monkey.filled.1.fa --repeat /panfs/GAG/
liubinghang/qc_assembly/crab_eating/500bac/test/repeat.list --step 12345 --outdir /panfs/GAG/liubinghang/qc_assembly/crab_e
ating/500bac/test/bac &

Please remember scaffold is one, and bac is two
=cut

use strict;
use Getopt::Long;

my $path="/ifs1/GAG/assemble/liubinghang/package/qc_assembly/BAC";

my ($bac_len,$scaf_len,$bac_seq,$scaf_seq);
my ($scaf_gap_file,$bac_gap_file,$bac_repeat_file,$scaf_repeat_file,$bac_gene_file,$scaf_gene_file);
my ($step);
my ($Verbose,$Help,$Outdir);
GetOptions(
        "len2:s"=>\$bac_len,
        "len1:s"=>\$scaf_len,
        "seq2:s"=>\$bac_seq,
        "seq1:s"=>\$scaf_seq,
        "step:i"=>\$step,
        "gap1:s"=>\$scaf_gap_file,
        "gap2:s"=>\$bac_gap_file,
        "repeat2:s"=>\$bac_repeat_file,
	"repeat:s"=>\$scaf_repeat_file,
        "gene2:s"=>\$bac_gene_file,
	"gene:s"=>\$scaf_gene_file,
        "outdir:s"=>\$Outdir,
        "verbose"=>\$Verbose,
        "help"=>\$Help
);
$step ||= 1;
$Outdir ||= ".";

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my (%h,%l,%b);

if($step =~ /1/){
	print "Step1 Preparing...\n";
	`fastaDeal.pl  --attr id:len $bac_seq >$Outdir/bac_seq.len`;
	`fastaDeal.pl  --attr id:len $scaf_seq >$Outdir/scaf_seq.len`;
	#`perl /panfs/GAG/liubinghang/bin/fasta_bin/caculate_gapsize.pl $bac_seq >$Outdir/bac.fa.contig.position`;
	#`perl /panfs/GAG/liubinghang/bin/fasta_bin/caculate_gapsize.pl $scaf_seq >$Outdir/scaf.fa.contig.position`;
	print "step1 is over!\n";
}

$bac_len="$Outdir/bac_seq.len";
$scaf_len="$Outdir/scaf_seq.len";
#$bac_gap_file="$Outdir/bac.fa.contig.position";
#$scaf_gap_file="$Outdir/scaf.fa.contig.position";

if($step =~ /2/){
	print "step 2: compare the genome and bac sequence!\n";
	open WR,">$Outdir/nucmer.sh";
	print WR "nucmer -p $scaf_seq $bac_seq $scaf_seq;";
	print WR "delta-filter -1 $scaf_seq.delta > $scaf_seq.delta.filter;";
	print WR "show-coords -q -T -H -d $scaf_seq.delta.filter > $scaf_seq.delta.filter.coors;";
	print WR "perl /nas/GAG_02/fanw/GACP/GACP-7.0/06.evolution_analysis/synteny_drawing_parallel/bin/convert_nucmer_to_synteny.pl $scaf_seq.delta.filter.coors > $Outdir/scaf_bac.delta.filter.coors.list;\n";
	close WR;
	print "Begin qsub, please waiting...\n";
	`nohup qsub-sge.pl --maxjob 1 --reqsub --jobprefix nm $Outdir/nucmer.sh`;
	`rm $scaf_seq.delta.filter.coors`;
	`rm $scaf_seq.delta.filter`;
	`rm $scaf_seq.cluster`;
	`rm $scaf_seq.delta`;
	print "step2 is over!\n";
}
##get the raw pair list
open (BL,$bac_len)||die "fail open $bac_len \n";
while(<BL>){
        my @t=split /\t/;
        $h{$t[0]}=1;
	$l{$t[0]}=1;
}
close BL;

my $list="$Outdir/scaf_bac.delta.filter.coors.list";
if($step =~ /3/){
	mkdir ("$Outdir/pair")  unless(-d "$Outdir/pair");
	foreach my $k (sort keys %h){
	        print "Finding bac id:  $k\n";
	        `more $list|grep '$k' >$Outdir/pair/$k.list`;
		`perl $path/get_all_pair.pl $Outdir/pair/$k.list $bac_len $scaf_len >$Outdir/pair/$k.log`;
	}
	`cat $Outdir/pair/*.log >$Outdir/pair.list`;
#	`perl /panfs/GAG/liubinghang/package/qc_assembly/BAC/BAC_coverage.pl $bac_len $Outdir/pair.list >$Outdir/total.coverage`;
	print "step 3 is over1\n";
}

my $pair_list="$Outdir/pair.list";
my %len;
open (PL,$pair_list) ||die "fail open $pair_list!\n";
while(<PL>){
        chomp;
	my @t=split /\t/;
        $b{$t[0]}=$t[3];
	$len{$t[0]}=$t[8];
}
close PL;

##get the pair seq, blast and solar.
if($step =~ /4/){
	print "Step 4: detail blast each pair,solar the result, get the statistic result\n";
        `perl $path/BAC_coverage.pl $bac_len $Outdir/pair.list >$Outdir/total.coverage`;
	mkdir "$Outdir/tmp" unless(-d "$Outdir/tmp");
	mkdir "$Outdir/output"  unless(-d "$Outdir/output");
	`perl $path/get_pair_seq.pl $pair_list $scaf_seq $bac_seq $Outdir >$Outdir/blast.sh`;
	`nohup qsub-sge.pl --maxjob 5 --reqsub --jobprefix nm $Outdir/blast.sh`;
	print "blast finish!\n";
	foreach my $k (sort keys %b){
        	`perl $path/indel/solar.pl $Outdir/tmp/$k/$k.$b{$k}.m8 1>$Outdir/tmp/$k/$k.$b{$k}.blast.solar 2>$Outdir/tmp/$k/$k.$b{$k}.blast.seq`;
		`perl $path/solar2list.pl  $Outdir/tmp/$k/$k.$b{$k}.blast.solar > $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list`;
		
		`perl $path/indel/extract_block.pl $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list $Outdir/tmp/$k/$k.$b{$k}.blast.seq >$Outdir/tmp/$k/$k.$b{$k}.blast.seq.block`;
		`perl $path/indel/get_gapregion2.pl $Outdir/tmp/$k/$k.$b{$k}.blast.seq.block`;
		#`perl $path/BAC_stat.pl $Outdir/bac_seq.len $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list >$Outdir/tmp/$k.stat`;
		#`perl $path/get_gap_region.pl $Outdir/tmp/$k/$b{$k}.fa $Outdir/tmp/$k/$k.fa $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.chr $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list >$Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.gap.block`;
#		`more $scaf_repeat_file | grep '$k' >$Outdir/tmp/scaf_repeat.list`;
#		`perl /panfs/GAG/liubinghang/package/qc_assembly/BAC/gapN_caculate2.pl $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.indel $Outdir/tmp/$k/$b{$k}.fa $Outdir/tmp/$k/$k.fa $Outdir/tmp/scaf_repeat.list >$Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.indel.full`;
		`perl $path/gapN_caculate.pl $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.indel $Outdir/tmp/$k/$b{$k}.fa $Outdir/tmp/$k/$k.fa >$Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.indel.full`;
		#`mv $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.indel.full $Outdir/tmp/`;
		`mv $Outdir/tmp/$k/$k.$b{$k}.blast.seq.block.indel $Outdir/output/`;
		`cp $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list $Outdir/tmp/`;	
	}
	`cat $Outdir/tmp/*.solar.list >$Outdir/output/all.solar.list`;
	`rm $Outdir/tmp/*.solar.list`;
	#`cat $Outdir/tmp/*.stat >$Outdir/output/total.statistic`;
	#`cat $Outdir/tmp/*.indel.full >$Outdir/output/total.indel.full`;
#	`perl /panfs/GAG/liubinghang/package/qc_assembly/BAC/get_table.pl $Outdir/output/total.statistic $Outdir/pair.list.coverage $bac_len >$Outdir/output/bac_coverage.table`;
	#`rm $Outdir/tmp/*.stat`;
	print "step4 is over!\n";
}

##draw the pictures in png

if($step =~ /5/){
	`sort $Outdir/output/all.solar.list +8 -g >$Outdir/output/all.solar.list.sort`;
	`perl $path/BAC_stat2.pl $Outdir/bac_seq.len $Outdir/output/all.solar.list.sort >$Outdir/output/coverage.table`;
	foreach my $k (sort keys %b){
		 print "Drawing scaffold: $k\n";
		`perl /panfs/GAG/liubinghang/qc_assembly/crab_eating/500bac/bin/bac.pl $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list >$Outdir/tmp/$k/$k.$b{$k}.log`;
	        `perl /panfs/GAG/liubinghang/bin/fasta_bin/caculate_gapsize.pl $Outdir/tmp/$k/$b{$k}.fa >$Outdir/tmp/$k/contig2.list`;
	        `perl /panfs/GAG/liubinghang/bin/fasta_bin/caculate_gapsize.pl $Outdir/tmp/$k/$k.fa >$Outdir/tmp/$k/contig1.list`;
		`perl $path/get_gap_region.pl $Outdir/tmp/$k/$b{$k}.fa $Outdir/tmp/$k/$k.fa $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.chr $Outdir/tmp/$k/contig1.list $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list >$Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.gap.block`;
		`mv $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.gap.block $Outdir/output/`;
		#`more $scaf_gap_file |grep '$k' >$Outdir/tmp/contig1.list`;
		#`more $bac_gap_file | grep '$b{$k}' >$Outdir/tmp/contig2.list`;
		`more $bac_repeat_file | grep '$b{$k}' >$Outdir/tmp/bac_repeat.list`;
		`more $scaf_repeat_file | grep '$k' >$Outdir/tmp/scaf_repeat.list`;
		`more $bac_gene_file |grep '$b{$k}' >$Outdir/tmp/gene.list`;
		`more $scaf_gene_file |grep '$k' >$Outdir/tmp/gene2.list`;
#		`perl /panfs/GAG/liubinghang/package/draw_svg/draw_scaf_bac2.pl --outdir $Outdir/tmp/$k --repeat2 $Outdir/tmp/scaf_repeat.list --gene2 $Outdir/tmp/gene2.list --gene $Outdir/tmp/gene.list --repeat $Outdir/tmp/bac_repeat.list --len $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.chr --list $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list --contig1 $Outdir/tmp/contig1.list --contig2 $Outdir/tmp/contig2.list --scale 10000 --resolution 100`;
		`perl /panfs/GAG/liubinghang/package/draw_svg/draw_scaf_bac2.pl --outdir $Outdir/tmp/$k --gene $Outdir/tmp/gene.list --repeat $Outdir/tmp/bac_repeat.list --len $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list.chr --list $Outdir/tmp/$k/$k.$b{$k}.blast.solar.list --contig1 $Outdir/tmp/$k/contig1.list --contig2 $Outdir/tmp/$k/contig2.list --scale 10000 --resolution 100`;
		`perl /panfs/GAG/liubinghang/package/draw_svg/svg2xxx.pl $Outdir/tmp/$k/$k\_$len{$k}.$b{$k}.parallel.svg -type png -memory 3G`;
		`mv $Outdir/tmp/$k/$k\_$len{$k}.$b{$k}.parallel.svg $Outdir/output/`;
		`mv $Outdir/tmp/$k/$k\_$len{$k}.$b{$k}.parallel.png $Outdir/output/`;
	}
	print "step 5 is over!\n";
}

