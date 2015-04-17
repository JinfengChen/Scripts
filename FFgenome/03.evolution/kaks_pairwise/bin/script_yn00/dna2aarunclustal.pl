#!/usr/bin/perl

my $directory=$ARGV[0];
my @file;
my $file;

###read in the directory and translate the dna file into aa file.

opendir DH, $directory or die "Can't open dir";
foreach $file (readdir DH) {
     $/=">";
     if ($file=~/(.*)\.fa/) {
                my $temp1="$directory/$file";
                my $temp2="$directory/$1aa.fa";
		open INFILE, "$directory/$file" or die "can not open my dna files";
		open OUTFILE, ">$directory/$1aa.fa" or die "can not open my aa file";
		while (my $read=<INFILE>) {
			  chomp $read;
                          my $name;
			  my $seq;
			  my $aaseq;
			  my @lines=split("\n",$read);
                          chomp @lines; 
			  $name=shift @lines;
			  print "$name\n";
			  $seq=join("",@lines);
			  if ($seq=~/^\s$/){
		              next;
			  }elsif ($seq=~/\w+/) {
                              if (dna2peptide($seq) eq "Bad"){
                                  `rm $temp1 $temp2`;
                                  next;
                              }else{
                                  $aaseq=dna2peptide($seq);
                                  print OUTFILE ">$name\n$aaseq\n";
                              }
			}
		}
		close OUTFILE;
		close INFILE;
	}
	
}
closedir DH;

##run clustalw to align these aa file.

opendir DIR, $directory or die "can not open dir1";
          foreach my $file1 (readdir DIR) {
                if ($file1=~/(.*aa)\.fa/) {
		    my @file1;
		    push (@file1, $file1);
		    system "clustalw -infile=$directory/$file1 -outorder=input";
                    system "rm $directory/$file1 $directory/$1.dnd";
                }
          }
closedir DIR;


system "rm $directory/*aa.fa";
system "rm $directory/*aa.dnd";

## sub to tranlate dna to peptide.
sub dna2peptide {

    my($dna) = @_;

    use strict;
    use warnings;
    

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        my $temp=codon2aa(substr($dna,$i,3)); 
        if ($temp eq "Bad"){
           return "Bad";
        }else{
           $protein .= codon2aa( substr($dna,$i,3) );
        }
    }
    #print "in sub, $protein\n";
    return $protein;
}


sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
        print STDERR "Bad codon \"$codon\"!!\n";
        return "Bad";    
    }
}



