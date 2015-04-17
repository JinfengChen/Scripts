#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

# assume short read fastq with 4 lines per record!

my $size = 1_000_000;
my $outdir;

if (!defined @ARGV){
  print "Please Provide 
-s int The number of sequences per file [1_000_000]
-o dir The name of the directory to output the files
followed by a list of files to be split

Usage:
./fastq_split -o split_fq ~/somedir/someRandom.fq
./fastq_split -o split_fq ~/somedir/*fq

";
exit;
}

GetOptions('s|size:i'  => \$size,
	   'o|outdir:s' => \$outdir,
	   );
my @files = @ARGV;

if (!-e $outdir){
  system ("mkdir -p $outdir");
}

for my $file ( @files ) {
    $file = File::Spec->rel2abs($file);
    my ($vol,$dir,$fname) = File::Spec->splitpath($file);
    my $odir = $outdir || $dir;

    my $fh;
    if( $fname =~ /(\S+)\.(gz|bz2)$/) {
	open($fh => "$uncompress{$2} $file |") ||die $!;
	$fname = $1;
    } else {
	open($fh => $file) || die "$file $!";
    }
    my @name = split(/\./,$fname);
    my ($ext) = pop @name;
    my $f = join(".",@name);    
    my $i = 0;
    #open(my $ofh => sprintf(">$odir/%s.p%02d.%s",$f,$i++,$ext)) || die $!;
    open(my $ofh => sprintf(">$odir/p%02d.%s.%s",$i++,$f,$ext)) || die $!;
    my $n = 0;
    while(<$fh>) {
	unless( /^@/ ) {
	    chomp;
	    die("out of register, got line:\n$_\n  but expected line to start with '\@'\n");
	}
	print $ofh $_;
	for ( 0..2) {
	    $_ = <$fh>;
	    print $ofh $_;
	}
	if( ++$n >= $size ) {
	    close($ofh);
	    #open($ofh => sprintf(">$odir/%s.p%02d.%s",$f,$i++,$ext)) || die $!;
    	    open($ofh => sprintf(">$odir/p%02d.%s.%s",$i++,$f,$ext)) || die $!;
	    $n= 0;
	}
    }
    close($ofh);
}


