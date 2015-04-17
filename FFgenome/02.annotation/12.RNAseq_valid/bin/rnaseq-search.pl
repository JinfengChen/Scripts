#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use File::Basename;
use Cwd 'abs_path';

pod2usage(0) unless $ARGV[0];
my ($_help, $overlap, $dir, $fasta, $ext, $border, $test, $format,
    $align, $program, $opt, $outdir);
GetOptions('h|help'   => \$_help,
	   'overlap'  => \$overlap,
	   'd|dir'      => \$dir,
	   'fasta'    => \$fasta,
	   'ext=s'    => \$ext,
	   'format=s' => \$format,
	   'border=s' => \$border,
	   'test=s'     => \$test,
	   'align'    => \$align,
	   'program=s' => \$program,
	   'opt=s'    => \$opt,
	   'outdir=s'    => \$outdir,
	   ) or die;
pod2usage(2) if $_help;

$format ||= 'eland';
my $targetfile = shift;
my $mapfile = shift;
my @mapfiles = @ARGV;

$ext ||= '_export.txt';
$border ||= 500;
if ($outdir) {$outdir = abs_path($outdir);}

if ($align) {
    opendir DIR, "$mapfile";
    foreach my $file (readdir(DIR)) {
	if ($file !~ /$ext$/) { next; }
	my $outfile = basename($mapfile);
	$file =~ s/$ext$//;
	$outfile = "$outfile-$file";
	my $job = create_job($outdir, $program, $targetfile, "$mapfile/$file$ext", $outfile);
	print STDERR "$outdir $job\n";
	system("qsub -V -l virtual_free=4G -o $job.err -e $job.out $job");
    }
    exit;
}

my ($gene_coords, $exon_coords, $intron_coords, $genes) 
    = get_coord_from_gff($targetfile);

if (!$dir) { 
    foreach my $file ($mapfile, @mapfiles) {
	my $outfile = basename($file);
	my ($count, $pos) = count_mappings($file); 
	print_gene_counts($count, "$outfile.gene");
	print_pos($pos, "$outfile.pos");
    }
}else {
    opendir DIR, "$mapfile";
    foreach my $file (readdir(DIR)) {
	if ($file !~ /$ext$/) { next; }
	my $outfile = basename($mapfile);
	$file =~ s/$ext$//;
	my ($count, $pos) = count_mappings("$mapfile/$file$ext");
	print_gene_counts($count, "$outfile-$file.gene");
	print_pos($pos, "$outfile-$file.pos");
    }
}
exit;

# create all job files
sub create_jobs {
    my ($out_dir, $files) = @_;
    my @jobs;
    for (my $i = 1; $i < @$files; $i++) {
	push @jobs, create_job($out_dir, $files->[$i]);
    }
    return \@jobs;
}

# create a job file
sub create_job {
    my ($out_dir, $program, $target, $query, $outfile, $param) = @_;
    my $file = basename($query);
    $outfile .= ".".basename($program);
#    print "$query\n"; exit;

    my $command = '';
    if ($program =~ /soap/) {
	$command .= "$program -D $target -a $query -o $out_dir/$outfile -r 0";
    }

    my $job = "$out_dir/job.$file.sh";
    open (BLAST, ">$job");

    print BLAST <<END;
#!/bin/sh

cd $ENV{PWD}
$command

END

    close BLAST; 
    system("chmod a+x $job");
#    print STDERR "chmod a+x $job\n";
    return $job;
}

sub count_mappings {
    my ($mapfile) = @_;
    print STDERR "# $mapfile\n"; 
    my %map_count;
    my %pos_on_target;
    open MAPPING, "<$mapfile" or die;
    
    my $mapped = 0;
    my $on_gene = 0;
    my $multi = 0;
    while (1) {
	my $query;
	if ($format eq 'soap') { $query = read_mapping_soap(\*MAPPING); }
	elsif ($format eq 'bowtie') {$query = read_mapping_bowtie(\*MAPPING); }
        elsif ($format eq 'bed') {$query = read_mapping_bed(\*MAPPING); }
	else {$query = read_mapping(\*MAPPING);}
	if (!$query) {last;}
	if ($test && scalar(keys %map_count) > $test) {last;}
	$mapped++;
	print STDERR join(' ', @$query), "\n";
	my $all_genes = $gene_coords->{$query->[1]};
	if (!$all_genes) {print STDERR " no genes\n";next;}
	my $index = bsearch($query->[3], $all_genes);
#    print join(' ', @{$all_genes->[$index-2]}), "\n";
#    print join(' ', @{$all_genes->[$index-1]}), "\n";
#    print join(' ', @{$all_genes->[$index]}), "\n";

	my @good_genes;
	my @opp_genes;
	my @border_genes;
	foreach my $i ($index..@{$all_genes}-1) {
	    my $g = $all_genes->[$i];
	    if ($g->[2] > $query->[4]+$border) {last;}
	    if ($g->[2]<=$query->[4] && $g->[3]>=$query->[3]) {
		if ($g->[1] == $query->[2]) {push @good_genes, $g; print STDERR ' same: ', join(' ', @$g), "\n";}
		else {push @opp_genes, $g;print STDERR ' opp: ', join(' ', @$g), "\n";}
	    }elsif ($g->[2]-$border<$query->[3] && $g->[3]+$border>$query->[4]) {
		push @border_genes, $g;print STDERR ' border: ', join(' ', @$g), "\n";
	    }
	}
	my $index1 = $index-5;
	if ($index1 < 0) {$index1 = 0;}
	foreach my $i ($index1..$index-1) {
	    my $g = $all_genes->[$i];
#	    if ($g->[2] > $query->[4]+$border) {last;}
	    if ($g->[2]<=$query->[4] && $g->[3]>=$query->[3]) {
		if ($g->[1] == $query->[2]) {push @good_genes, $g;print STDERR ' same: ', join(' ', @$g), "\n";}
		else {push @opp_genes, $g;print STDERR ' opp: ', join(' ', @$g), "\n";}	    
	    }elsif ($g->[2]-$border<$query->[3] && $g->[3]+$border>$query->[4]) {
		push @border_genes, $g;print STDERR ' border: ', join(' ', @$g), "\n";
	    }
	}
	
	if (@good_genes || @opp_genes) {$on_gene++;}
	elsif (@border_genes) {}
	else {print STDERR " no genes\n"; next}
	if (@good_genes && @opp_genes) {print STDERR " opposite genes\n";}
	my $on_exon = 0;
	foreach my $g (@good_genes, @opp_genes) {
	    my $g_id = $g->[0];
#	    print " $g_id exons:\n";
	    my $in_exon = 0;
	    my $over_exon = 0;
	    foreach my $exon (@{$exon_coords->{$g_id}}) {
#		print "  exon: ", join(' ', @$exon), "\n";
		if ($query->[3]>= $exon->[2] && $query->[4]<=$exon->[3]) {
		    $in_exon = 1;
		    $map_count{$g_id}{in_exon}++;
		    $map_count{$g_id}{exons}{$exon}=$exon;
		    print STDERR "  in exon", ($on_exon?' multi':''), "\n"; 
		    last;
		}elsif ($query->[3]<$exon->[2] && $query->[4]>=$exon->[2]
			||$query->[3]<=$exon->[3] && $query->[4]>$exon->[3]) {
		    $over_exon = 1;
		}
	    }
	    if ($over_exon == 1 && $in_exon == 0) {
		$map_count{$g_id}{over_exon}++;
		print STDERR "  over exon", ($on_exon?' multi':''), "\n"; 
	    }
	    if (($in_exon || $over_exon) && $on_exon) {$multi++;}
	    if ($in_exon || $over_exon) {
		$on_exon=1;
		check_pos($query, $g, \%pos_on_target);
	    }
#	    else {print "  NOT on exon\n";}
	}
	if ($on_exon) {next;}
	
	my $opp_exon = 0;
	foreach my $g (()) {
	    my $g_id = $g->[0];
	    foreach my $exon (@{$exon_coords->{$g_id}}) {
		if ($query->[3]<= $exon->[3] && $query->[4]>=$exon->[2]) {
		    $map_count{$g_id}{opp_exon}++;
		    $opp_exon = 1;
		}
	    }
	}
#    if ($opp_exon) {next;}
	my $in_intron = 0;
	foreach my $g (@good_genes, @opp_genes) {
	    my $g_id = $g->[0];
#	    print " $g_id introns:\n";
#	foreach my $intron (@{$intron_coords->{$g_id}}) {
#	    print "  intron: ", join(' ', @$intron), "\n";
#	    if ($query->[3]>=$intron->[2] && $query->[4]<=$intron->[3]) {
	    $map_count{$g_id}{in_intron}++;
	    if ($in_intron) {$multi++;}
	    print STDERR "  in intron", ($in_intron?' multi':''), "\n";
	    $in_intron=1;
#		last;
#	    }
#	}
	}
	if ($in_intron) {next;}

	my $on_border = 0;
	foreach my $g (@border_genes) {
	    my $g_id = $g->[0];
#	    print " $g_id borders:\n";
	    $map_count{$g_id}{on_border}++;
	    print STDERR "  on border", ($on_border?' multi':''), "\n";
	    $on_border = 1;
	    check_pos($query, $g, \%pos_on_target);
	}
#	print "on gene: $on_exon, $in_intron\n";
    }
    close MAPPING;
    
    print STDERR "==mapped: $mapped, on gene: $on_gene, extra mappings: $multi\n";

    return (\%map_count, \%pos_on_target);
}

sub print_gene_counts {
    my ($map_count, $outfile) = @_;
    open OUT, ">$outfile" or die;
    my %map_count = %$map_count;
    my @map_counts = map {
	my $num=0;
	foreach my $k (keys %{$map_count{$_}}) {
	    if ($k ne 'in_exon' && $k ne 'over_exon') {next;}
	    $num += $map_count{$_}{$k}||0;
	}
	[$_, 
	 $map_count{$_}{in_exon}||0,
	 $map_count{$_}{over_exon}||0,
	 $map_count{$_}{opp_exon}||0,
	 $map_count{$_}{in_intron}||0,
	 $genes->{$_},
	 $num,
	 $map_count{$_}{on_border}||0,
	 $map_count{$_}{exons}
	 ]
	 } keys %map_count;
    @map_counts = sort {$b->[6] <=> $a->[6]
			    or $b->[1] <=> $a->[1]
			    or $b->[2] <=> $a->[2]
			    or $b->[7] <=> $a->[7]
			} @map_counts;
    my $num = 0;
    my $b =0;
    print OUT join("\t", 'gene_id', 'in_exon', 'over_exon', 'opp_exon', 'in_intron', 'gene_name', 'on_exon_per_kb', "on_border_$border"), "\n";
    foreach my $g (@map_counts) {
	my $exon_len = 0;
#	print join(' ', keys %{$g->[8]}), ', ', join(' ', values %{$g->[8]}), "\n";
	foreach my $e (values %{$g->[8]}) {
	    $exon_len += $e->[4];
	}
	print STDERR $g->[0], ' ', $g->[1], " exon len: ", $exon_len,"\n";
	$g->[6] = int($g->[6]/(($exon_len/1000)||1));
	$num += $g->[6];
	$b += $g->[7];
	if ($g->[6] == 0 && $g->[7]<2) {next;}
	print OUT join("\t", @$g[0..(@$g-2)]), "\n";
    }
    close OUT;
    print STDERR "==all mapped to gene: $num, border: $b\n===\n";
}

sub print_pos {
    my ($p_on_t, $outfile) = @_;
    my %pos5;
    my %pos3;
    my @pos = (0)x101;
    foreach my $g (keys %$p_on_t) {
	print STDERR "gene: $g\n";
	my $p_on_g = $p_on_t->{$g};
	my $g_len = $p_on_g->[0];
	my $middle = int($g_len/2);
	foreach my $p (keys %{$p_on_g->[1]}) {
	    if ($p>$g_len) {die "pos err: $p, $g_len\n";}
	    if ($p <= 300 && $p <= $middle) {
		$pos[$p] += $p_on_g->[1]->{$p};
	    }elsif ($g_len-$p < 300 && $p > $middle) {
		$pos[$p+1300] += $p_on_g->[1]->{$p};
	    }else {
		$pos[int(1000*($p-300)/($g_len-600)+0.5)+300] += $p_on_g->[1]->{$p};
	    }
	}
	foreach my $p (keys %{$p_on_g->[2]}) {
	    $pos5{$p} += $p_on_g->[2]->{$p};
	}
	foreach my $p (keys %{$p_on_g->[3]}) {
	    $pos3{$p+100} += $p_on_g->[3]->{$p};
	}
    }
    open OUT, ">$outfile" or die;
    foreach my $p (sort {$b<=>$a} keys %pos5) {
	print OUT "-$p\t", $pos5{$p}, "\n";
    }
    foreach my $i (0..@pos-1) {
	print OUT $i, "\t", $pos[$i], "\n";
    }
    foreach my $p (sort {$a<=>$b} keys %pos3) {
	print OUT $p, "\t", $pos3{$p}, "\n";
    }
    close OUT;
}

sub read_mapping {
    my $fh = shift;
    while (my $line = <$fh>) {
	chomp $line;
	my @items = split "\t", $line;
	my $seq = $items[8];
	my $qs = $items[9];
	my $chr = $items[11];
	if (!$chr) {next;}
	my $start = $items[12];
	my $end = $start+length($seq)-1;
	my $strand = $items[13] eq 'F' ? 1 : -1;
	my $align = $items[14];
	my $align_score = $items[15];
	my $filter = $items[21];
	return [$seq,$chr,$strand,$start,$end, $align, $align_score, $filter];
    }
}
sub read_mapping_bed {
    my $fh = shift;
    #print "Reading bed align ...\n";
    while (my $line = <$fh>) {
        chomp $line;
        my @items = split "\t", $line;
        my $seq = $items[3];
        #my $qs = $items[3];
        my $chr = $items[0];
        if (!$chr) {next;}
        my $start = $items[1];
        my $end = $items[2];
        my $strand = $items[5] eq '+' ? 1 : -1;
        my $align = "";
        my $align_score = "";
        my $filter = "Y";
        return [$seq,$chr,$strand,$start,$end, $align, $align_score, $filter];
    }
}

sub read_mapping_bowtie {
    my $fh = shift;
    while (my $line = <$fh>) {
	chomp $line;
	my @items = split "\t", $line;
#	print STDERR join(', ', @items), "\n";
	my $seq = $items[4];
	my $qs = $items[5];
	$items[2] =~ /(\S+)$/;
	my $chr = $1;
	if (!$chr) {next;}
	my $start = $items[3];
	my $end = $start+length($seq)-1;
	my $strand = $items[1] eq '+' ? 1 : -1;
	my $align = $items[6];
	my $align_score = 100;
	my $filter = 'Y';
	return [$seq,$chr,$strand,$start,$end, $align, $align_score, $filter];
    }
}

sub read_mapping_soap {
    my $fh = shift;
    while (my $line = <$fh>) {
	chomp $line;
	my @items = split "\t", $line;
	my $seq = $items[1];
	my $qs = $items[2];
	my $chr = $items[7];
	if (!$chr) {next;}
	if ($chr =~ /^\w+:\w+:(\w+)/) {$chr = $1;}
	my $start = $items[8];
	my $end = $start+length($seq)-1;
	my $strand = $items[6] eq '+' ? 1 : -1;
	my $align = $items[-1];
	my $align_score = $items[9];
	my $filter = 'Y';
	return [$seq,$chr,$strand,$start,$end, $align, $align_score, $filter];
    }
}

sub get_coord_from_fasta {
    my $seqin  = Bio::SeqIO->new(-file => "$targetfile" , '-format' => 'Fasta');
    my %gene_coords;
    my @locus_coords;
    for (my $i=0; my $seq = $seqin->next_seq(); $i++) {
	my ($gene_id) = split ':', $seq->id;
	print STDERR $gene_id, "\n";
#	print $seq->description, "\n";
	my %items = split /;\s*|=/, $seq->description;
	print $items{loc}, "\n"; 
#	push @locus_coords, 
    }
    
}

sub get_coord_from_gff {
    my $file = shift;
    my %gene_coords;
    my %trans;
    my %genes;
    my %exon_lens;
    my %description;
    my %exon_coords;
    my %intron_coords;
    my $exon_num=0;
    my $intron_num=0;
    open IN, "<$file" or die "can't open $file: $!";
    while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^\#\#\#$/) {last;}
	if ($line =~ /^\#\#/) {next;}
	my @items = split "\t", $line;
	my $chr = $items[0];
	my $type = $items[2];
	if ($type !~ /gene|mRNA|CDS|exon|intron/) {next;}
	my $start = $items[3];
	my $end = $items[4];
	my $strand = $items[6] eq '+' ? 1 : -1;
	my %desc = split /=|;/, $items[8];
	my $id = $desc{ID};
#	if ($type eq 'gene') {print "$id\n";}
	
	if ($type eq 'gene') {
	    my $name = $desc{Name} || $desc{'ID'};
	    $genes{$id} = $name;
#	    $description{$id} = 
	    push @{$gene_coords{$chr}}, [$id,$strand,$start,$end, $chr, $name];
	}elsif ($type =~ /mRNA|pseudogene|transcript/) {
	    my $gene_id = $desc{Parent};
	    $trans{$id} = $gene_id;
	}elsif ($type eq 'CDS') {
	    my $tran_id = $desc{Parent};
	    ($tran_id) = split ',', $tran_id;
	    push @{$exon_coords{$tran_id}}, [$chr, $strand, $start, $end, $end-$start+1];
	    $exon_num++;
	}elsif ($type eq 'intron') {
	    my $tran_id = $desc{Parent};
	    ($tran_id) = split ',', $tran_id;
	    push @{$intron_coords{$tran_id}}, [$chr,$strand,$start,$end];
	    $intron_num++;
	}
    }
    close IN;
    my $num = 0;
    foreach my $k (keys %gene_coords) {
	$gene_coords{$k} = [sort {$a->[2]<=>$b->[2]} @{$gene_coords{$k}}];
	$num += @{$gene_coords{$k}};
    }
    my %exons;
    foreach my $k (keys %exon_coords) {
#	print STDERR "$k\n";
	push @{$exons{$trans{$k}}}, @{$exon_coords{$k}};
#	$exon_lens{$trans{$k}} ||= 0;
	foreach my $e (@{$exon_coords{$k}}) {
	    $exon_lens{$trans{$k}} += $e->[4];
	}
    }
    my %introns;
    foreach my $k (keys %intron_coords) {
#	print STDERR "$k\n";
	push @{$introns{$trans{$k}}}, @{$intron_coords{$k}};
    }
    print STDERR '# ', scalar(keys %gene_coords), " chrs, $num genes\n";
    print STDERR "# $exon_num exons, $intron_num introns\n";
    return (\%gene_coords, \%exons, \%introns, \%genes);
}


# positions on gene: 1 -- gene_len
# positions before gene: ... -- -1 (use positive number here)
# positions after gene: 1 -- ... 
sub check_pos {
    my ($query, $tg, $p_on_t) = @_;
    $p_on_t->{$tg->[0]} ||= [];
    my $tpos = $p_on_t->{$tg->[0]};
    $tpos->[0] = $tg->[3]-$tg->[2]+1;
    for (my $i = $query->[3]; $i <= $query->[4]; $i++) {
	if ($i < $tg->[2]) {
	    $tpos->[2]->{$tg->[2]-$i}++;
	}elsif ($i > $tg->[3]) {
	    $tpos->[3]->{$i-$tg->[3]}++;
	}else {
	    $tpos->[1]->{$i-$tg->[2]+1}++;
	}
    }
}

sub get_ids {
    my $idfile = shift;
    open(ID, "<$idfile");
    my @ids = <ID>;
    close ID;
    my $ids = join '', @ids;
    $ids =~ s/^\s+//;
    $ids =~ s/\s+$//;
    @ids = split /\s+/, $ids;
    close ID;
    return \@ids;
}

sub get_keywords {
    my $keyword = shift;
    open KEY, "<$keyword";
    my @keys;
    while (my $key = <KEY>) {
	$key =~ s/^\s+//;
	$key =~ s/\s+$//;
	push @keys, $key;
    }
    close KEY;
    return \@keys;
}

# search array of gene position a for given position x (integer)
# return index where found or the one on right side if not found
sub bsearch {
    my ($x, $a) = @_;            # search for x in array a
    my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
    my $i;                       # index of probe
    while ($l <= $u) {
	$i = int(($l + $u)/2);
#	print($i, "\n");
#	print join(' ', $a->[$i]), "\n";
	if ($a->[$i]->[2] < $x) {
	    $l = $i+1;
	}
	elsif ($a->[$i]->[2] > $x) {
	    $u = $i-1;
	} 
	else {
	    return $i; # found
	}
    }
    return $l;         # not found
}
sub bsearch1 {
    my ($x, $a) = @_;            # search for x in array a
    my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
    my $i;                       # index of probe
    while ($l <= $u) {
	$i = int(($l + $u)/2);
#	print($i, "\n");
#	print $a->[$i];
	if ($a->[$i] lt $x) {
	    $l = $i+1;
	}
	elsif ($a->[$i] gt $x) {
	    $u = $i-1;
	} 
	else {
	    return $i; # found
	}
    }
    return -1;         # not found
}

# ----------------------------------------------------

=head1 NAME

fasta-subset.pl

=head1 SYNOPSIS

    rnaseq-search target-file mapping-file -help -overlap

Options:

  -h|--help        Show brief help and exit

=head1 DESCRIPTION


=head1 SEE ALSO

perl.

=head1 AUTHOR

Chengzhi Liang E<lt>liang@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

