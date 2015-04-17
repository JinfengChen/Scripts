#!/usr/bin/perl


=head1 NAME
	
protein_map_genome.pl  --  the pipeline of mapping protein to genome for draft sequence.

=head1 DESCRIPTION

This script invokes tblastn, genBlastA, and genewise, and finally
give result in gff3 format.

The total time sometimes depends on the some tasks of genewise,
which becomes very slowly when the protein and DNA sequence is
too large. 

Alert that in this program, "|"  in protein name is not allowed.

=head1 Version

  Author: Huang Quanfei,huangqf@genomics.org.ch
  Version: 3.0  Date: 2009-7-10

=head1 Usage

    perl protein_map_genome.pl [options] protein.fa genome.fa
    --cpu <int>	         set the cpu number to use in parallel, default=3   
    --run <str>          set the parallel type, qsub, or multi, default=qsub
    --outdir <str>       set the result directory, default="./"
    --tophit <num>       select best hit for each protein, default no limit
    --blast_eval <num>   set eval for blast alignment, default 1e-5
    --align_rate <num>   set the aligned rate for solar result, default 0.01
    --extend_len <num>   set the extend length for genewise DNA fragment
    --step <num>         set steps for running, default 1. Set 1234 to run all steps, 1blast, 2solar,3genewise,4gff.
    --queue <str>        set the queue ,default no
    --lines <int>        set the --lines option of qsub-sge.pl,required.
    --verbose            output verbose information to screen  
    --help               output help information to screen  

=head1 Example

  perl ../bin/protein_map_genome.pl -cpu 100 -verbose -lines 100 -step 1234 ../input/Drosophila_melanogaster_protein.1000.fa ../input/test_chr_123.seq &


=cut


use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$Bin/../lib";
use Data::Dumper;
use File::Basename qw(basename dirname); 
use YAML qw(Load Dump);
use List::Util qw(reduce max min);
use Collect qw(Read_fasta);

##get options from command line into variables and set default values
my ($Cpu,$Run,$Outdir);
my ($Blast_eval,$Align_rate,$Extend_len,$Step);
my ($Cpu,$Verbose,$Help);
my $Queue;
my $Line;
my $Tophit;

GetOptions(
	"lines:i"=>\$Line,
	"cpu:i"=>\$Cpu,
	"run:s"=>\$Run,
	"outdir:s"=>\$Outdir,
	"blast_eval:s"=>\$Blast_eval,
	"align_rate:f"=>\$Align_rate,
	"extend_len:i"=>\$Extend_len,
	"step:s"=>\$Step,
	"cpu:i"=>\$Cpu,
	"queue:s"=>\$Queue,
	"tophit:i"=>\$Tophit,
	"verbose!"=>\$Verbose,
	"help!"=>\$Help
);
$Blast_eval ||= 1e-5;
$Align_rate ||= 0.25;
$Extend_len ||= 500;
$Step ||= '1';
$Cpu ||= 3;
$Run ||= "qsub";
$Outdir ||= ".";
my $Queue_para=(defined $Queue)?"-q $Queue":'';
#die `pod2text $0` if (@ARGV == 0 || $Help || (not defined $Line ));
die `pod2text $0` if (@ARGV == 0 || $Help);
my $Qr_file = shift;
my $Db_file = shift;

my %Pep_len;
read_fasta($Qr_file,\%Pep_len);

my %Chr_len;
read_fasta($Db_file,\%Chr_len);

$Outdir =~ s/\/$//;
mkdir($Outdir) unless(-d $Outdir);

my $Qr_file_basename = basename($Qr_file);
my $genewise_dir = "$Outdir/$Qr_file_basename.genewise";

my $tblastn_shell_file = "$Outdir/$Qr_file_basename.tblastn.shell";
my $genewise_shell_file = "$Outdir/$Qr_file_basename.genewise.shell";

my @subfiles;

##use YAML format to set parameters for blastall, solar, filter, and genewise programs
my $Param = Load(<<END);
blastall:
  -p: tblastn
  -e: $Blast_eval
  -F: F
  -m: 8
filter-solar:
  score: 25
  align_rate: $Align_rate
  extent: $Extend_len
genewise:
  -genesf:
  -gff:
  -sum:
END
print STDERR Dump($Param) if($Verbose);

if ($Step =~/1/){
	##format the database for tblastn
	print STDERR  "\n\n$Bin/blast/formatdb -i $Db_file -p F -o T\n"  if($Verbose);
	`$Bin/blast/formatdb -i $Db_file -p F -o T` unless (-f $Db_file.".nhr");

	##cut query file into small subfiles
	`perl $Bin/fastaDeal.pl -cutf $Cpu $Qr_file -outdir $Outdir`;
	@subfiles = glob("$Outdir/$Qr_file_basename.cut/*.*");

	##creat the tblastn shell file
	my $opt_blastall = join(" ",%{$Param->{blastall}});
	open OUT,">$tblastn_shell_file" || die "fail $tblastn_shell_file";
	foreach my $qrfile (@subfiles) {
		print OUT "$Bin/blast/blastall $opt_blastall -d $Db_file -i $qrfile -o $qrfile.blast; \n";
	}
	close OUT;

	print STDERR "run the tblastn shell file"  if($Verbose);
	if ($Run eq "qsub") {
		`perl $Bin/qsub-pbs.pl $Queue_para --maxjob $Cpu $tblastn_shell_file`;
	}
	if ($Run eq "multi") {
		`perl $Bin/multi-process.pl -cpu $Cpu $tblastn_shell_file`;
	}

	##cat together the tblastn result
	`cat $Outdir/$Qr_file_basename.cut/*.blast > $Outdir/$Qr_file_basename.blast`;
}

if ($Step =~/2/){
	print  STDERR "Run solar to conjoin HSPs and filter bad HSPs and redundance.\n"  if($Verbose);
	`perl $Bin/solar/solar.pl -a prot2genome2 -z -f m8 $Outdir/$Qr_file_basename.blast >$Outdir/$Qr_file_basename.blast.solar`;
	filter_solar("$Outdir/$Qr_file_basename.blast.solar");
	solar_to_table("$Outdir/$Qr_file_basename.blast.solar.filter","$Outdir/$Qr_file_basename.blast.solar.filter.table");
	`perl $Bin/genomic_cluster.pl  -overlap_percent 0.5 $Outdir/$Qr_file_basename.blast.solar.filter.table > $Outdir/$Qr_file_basename.blast.solar.filter.table.nonredundance`;
	`perl $Bin/fishInWinter.pl -bf table -ff table $Outdir/$Qr_file_basename.blast.solar.filter.table.nonredundance $Outdir/$Qr_file_basename.blast.solar.filter > $Outdir/$Qr_file_basename.blast.solar.filter.nr`;
	`rm -r $Outdir/$Qr_file_basename.blast.solar.filter.table*`;
}

if ($Step =~/3/){
	print "preparing genewise input directories and files\n" if ($Verbose);
	&prepare_genewise("$Outdir/$Qr_file_basename.blast.solar.filter.nr");

	print "run the genewise shell file\n" if ($Verbose);
	print  STDERR "running genewise\n"  if($Verbose);
	if ($Run eq "qsub") {
		`perl $Bin/qsub-pbs.pl $Queue_para  --maxjob $Cpu --lines $Line   $genewise_shell_file`;
	}
	if ($Run eq "multi") {
		`perl $Bin/multi-process.pl -cpu $Cpu $genewise_shell_file`;
	}
}

if ($Step =~/4/){
	print  STDERR "convert result to gff3 format\n"  if($Verbose);
	`for i in $Outdir/$Qr_file_basename.genewise/* ;do for j in \$i/*.genewise ;do cat \$j;done ;done >$Outdir/$Qr_file_basename.solar.genewise`;
	`perl $Bin/fastaDeal.pl -attr id:len $Qr_file >$Outdir/$Qr_file_basename.len`;
	`perl  $Bin/gw2gff.pl $Outdir/$Qr_file_basename.solar.genewise $Outdir/$Qr_file_basename.len >$Outdir/$Qr_file_basename.solar.genewise.gff`;
}

print  STDERR "All tasks finished\n"  if($Verbose);


##########################################################
################### Sub Routines ###################
#########################################################

##OsB000025-PA    476     1       476     -       Chr07frag1M     1000000 154122  157515  8       924     1,149;150,184;182,205;205,23
sub solar_to_table{
	my $file = shift;
	
	my $output;
	open IN, $file || die "fail";
	while (<IN>) {
		chomp;
		my @t = split /\t/;
		my $len = $t[3]-$t[2]+1;
		$output .= "$t[0]\t$t[5]\t$t[4]\t$t[7]\t$t[8]\t$len\n";
	}
	close IN;

	open OUT, ">$file.table" || die "fail";
	print OUT $output;
	close OUT;
}

##filter solar result, get parameters from globle $Param 
##################################################
sub filter_solar {
	my $infile = shift;
	my %solardata;
	my $output;
	
	open IN, "$infile" || die "fail $infile";
	while (<IN>) {
		chomp;
		s/^\s+//;
		my @t = split /\s+/;
		my $query = $t[0];
		my $score = $t[10];
		next if($score < $Param->{'filter-solar'}{score});
		my $query_size = $t[1];
		my $align_size;
		while ($t[11]=~/(\d+),(\d+);/g) {
			$align_size += abs($2 - $1) + 1;
		}
		next if($align_size / $query_size < $Param->{'filter-solar'}{align_rate});
	
		push @{$solardata{$query}},[$score,$_]; ## hits that better than cutoff
		
	}
	
	open OUT, ">$infile.filter" || die "fail $infile.filter";
	foreach my $query (sort keys %solardata) {
		my $pp = $solardata{$query};
		@$pp = sort {$b->[0] <=> $a->[0]} @$pp;
		for (my $i=0; $i<@$pp; $i++) {
			last if(defined $Tophit && $i>=$Tophit);
			my $query_Dup = "$query-D".($i+1);
			$pp->[$i][1] =~ s/$query/$query_Dup/ if ($i>0);
			print OUT $pp->[$i][1],"\n";
		}
	}
	close OUT;
	
}



##read sequences in fasta format and calculate length of these sequences.
sub read_fasta{
	my ($file,$p)=@_;
	open IN,$file or die "Fail $file:$!";
	$/=">";<IN>;$/="\n";
	while(<IN>){
		my ($id,$seq);
		#if ( /\S\s+\S/ ) {
		#	die "No descriptions allowed after the access number in header line of fasta file:$file!\n";
		#}
	#	if ( /\|/ ){
	#		die "No '|' allowed in the access number of fasta file:$file!\n";
	#	}
		
		if (/^(\S+)/){
			$id=$1;
		}else{
			die "No access number found in header line of fasta file:$file!\n";
		}
		if ( $id=~/\|/ ) {
			die "No '|' allowed in the access number of fasta file:$file!\n";
		}
		$/=">";
		$seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$p->{$id}=length($seq);
		$/="\n";
	}
	close IN;
}

##prepare data for genewise and make the qsub shell
####################################################


sub prepare_genewise{
	my $solar_file = shift;
	my @corr;

	open IN, "$solar_file" || die "fail $solar_file";
	while (<IN>) {
		s/^\s+//;
		my @t = split /\s+/;
		my $query = $t[0];
		my $strand = $t[4];
		my ($query_start,$query_end) = ($t[2] < $t[3]) ? ($t[2] , $t[3]) : ($t[3] , $t[2]);
		my $subject = $t[5];
		my ($subject_start,$subject_end) = ($t[7] < $t[8]) ? ($t[7] , $t[8]) : ($t[8] , $t[7]);
		push @corr, [$query,$subject,$query_start,$query_end,$subject_start,$subject_end,"","",$strand]; ## "query_seq" "subject_fragment"	
	}
	close IN;
	my %fasta;
	&Read_fasta($Qr_file,\%fasta);
	foreach my $p (@corr) {
		my $query_id = $p->[0];
		$query_id =~ s/-D\d+$//;
		if (exists $fasta{$query_id}) {
			$p->[6] = $fasta{$query_id}{seq};
		}
	}
	undef %fasta;
	my %fasta;
	&Read_fasta($Db_file,\%fasta);
	foreach my $p (@corr) {
		if (exists $fasta{$p->[1]}) {
			my $seq = $fasta{$p->[1]}{seq};
			my $len = $fasta{$p->[1]}{len};
			$p->[4] -= $Param->{'filter-solar'}{extent};
			$p->[4] = 1 if($p->[4] < 1);
			$p->[5] += $Param->{'filter-solar'}{extent};
			$p->[5] = $len if($p->[5] > $len);
			$p->[7] = substr($seq,$p->[4] - 1, $p->[5] - $p->[4] + 1); 
		}
	}
	undef %fasta;
	mkdir "$genewise_dir" unless (-d "$genewise_dir");
	my $subdir = "000";
	my $loop = 0;
	my $cmd;
	my $opt_genewise = join(" ",%{$Param->{genewise}});
	foreach my $p (@corr) {
		if($loop % 200 == 0){
			$subdir++;
			mkdir("$genewise_dir/$subdir");
		}
		
		my $qr_file = "$genewise_dir/$subdir/$p->[0].fa";
		my $db_file = "$genewise_dir/$subdir/$p->[0]_$p->[1]_$p->[4]_$p->[5].fa";
		my $rs_file = "$genewise_dir/$subdir/$p->[0]_$p->[1]_$p->[4]_$p->[5].genewise";
		
		open OUT, ">$qr_file" || die "fail creat $qr_file";
		print OUT ">$p->[0]\n$p->[6]\n";
		close OUT;
		open OUT, ">$db_file" || die "fail creat $db_file";
		print OUT ">$p->[1]_$p->[4]_$p->[5]\n$p->[7]\n";
		close OUT;

		my $choose_strand = ($p->[8] eq '+') ? "-tfor" : "-trev";
		#$cmd .= "$Bin/genewise $choose_strand $opt_genewise $qr_file $db_file > $rs_file 2> /dev/null;\n";
		$cmd .= "$Bin/genewise $choose_strand $opt_genewise $qr_file $db_file > $rs_file 2> /dev/null;\n";
                $loop++;
	}
	undef @corr;

	open OUT, ">$genewise_shell_file" || die "fail creat $genewise_shell_file";
	print OUT $cmd;
	close OUT;

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
	return ($pure_size);
}

