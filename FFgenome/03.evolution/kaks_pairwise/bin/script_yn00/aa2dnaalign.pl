#!/usr/bin/perl

my $directory=$ARGV[0];
my @aaalign;
my @dnaseq;


###aa align file and dna fasta file
opendir DIR1, $directory or die "can not open dir1";
          foreach my $file2 (readdir DIR1) {
                if ($file2=~/.*aa\.aln/) {
			push (@aaalign, $file2);
			#print "@aaalign\n";
                }elsif ($file2=~/.*\.fa/){
			push (@dnaseq, $file2);
			#print "@dnaseq\n";
	        }
          }
closedir DIR1;

@aaalign=sort @aaalign;
@dnaseq=sort @dnaseq;


## convert the aa align file to dna align file by using sub aa2dna.
for (my $i=0;$i<=$#dnaseq;$i++) {
    my $fileaa;
    my $filedna;
    my $description;
    my $fileout;
    $fileaa="$directory/$aaalign[$i]"; 
    $filedna="$directory/$dnaseq[$i]";
    $dnaseq[$i]=~/(.*)\.fa/;
    $description=$1;
    $fileout="$directory/$1.fas";
    #print "$fileaa\n$filedna\n$description\n$fileout\n";
    aa2dna($fileaa,$filedna,$description,$fileout);
    system ("rm $fileaa $filedna");
}
#system ("rm *.aln");

## sub to convert aa align to dna align.
sub aa2dna {
#chomp @_;
my ($file_aa,$file_dna,$description,$file_out)=@_;
my $i=0;
my $magic=0;
my $numofseq;
my $gene;
my @name;
my $seqnumber;
open INFILE, "$file_aa" or die "can not open my alignment file";
while (my $line=<INFILE>) {
     chomp $line;
     $i++;		   
     if ($i>=4 && $line=~/^\w/ && $magic==0) {###change 4 to 5 if the alignment is edit by genedoc generate a aln file.
         $numofseq++;
         my @word=split (" ", $line);
	 push (@name, $word[0]);
	 chomp @name;
     }elsif($i>=4){#####change 4 to 5 if the alignment is edit by genedoc generate a aln file.
         $magic=1;
     }		   
} 
close INFILE;
$i=0;

$/=">";
my %fasta;
my @keys;
my @values;
open INFILE0, "$file_dna" or die "can not open my fasta file";
     while (my $line=<INFILE0>) {
		   chomp $line;
		   my ($head,$sequence)=split (" ",$line);
		   chomp $head;
		   #print "$head\n";
		   chomp $sequence;
		   $fasta{$head}=$sequence; 
     }
     @keys=keys(%fasta);
     @values=values(%fasta);
     chomp @keys; 
     chomp @values;
close INFILE0;
$/="\n";
		
open INFILE1, "$file_aa" or die "can not open my alignment file";
open OUTFILE, ">$file_out" or die "can not open my outfile";              
while (my $line=<INFILE1>) {
      chomp $line;
      $i++;
      if ($i>=4 && $line=~/^\w/) {#####change 4 to 5 if the alignment is edit by genedoc generate a aln file.
          $seqnumber++;
          my @word=split (" ", $line);
	  my @aa=split("", $word[1]);
	  chomp @aa;
          my $len=@aa;
          push (@$seqnumber,@aa); 
      }elsif ($i>=4){#####change 4 to 5 if the alignment is edit by genedoc generate a aln file.
	  $seqnumber=0;
      }
}				
my @dnasequence;
for (my $j=1;$j<=$numofseq;$j++) {
     chomp @$j;
     print OUTFILE ">$name[$j-1]\n";
     @dnasequence=split("",$fasta{$name[$j-1]});
     my $length=@$j;
     my $codenumber;
     for (my $y=0;$y<$length;$y++) {### $y is the number of aa number in alignment including --.
             if ($$j[$y]=~/-/) {
		print OUTFILE "---";
             }elsif($$j[$y]=~/\w/){
	        for (my $z=0;$z<3;$z++) {
		    print OUTFILE "$dnasequence[3*$codenumber+$z]";
	        }
	        $codenumber++  ##count for codenumber
	     }
     }
     print OUTFILE "\n";
}
close OUTFILE;
close INFILE1;
for (my $j=1;$j<=$numofseq;$j++) {
    undef @$j;
}
}

