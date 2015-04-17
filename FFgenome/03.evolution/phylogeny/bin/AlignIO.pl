#!/usr/bin/perl
use Getopt::Long;
use Bio::AlignIO;

GetOptions (\%opt,"informat:s","outformat:s","inaln:s","outaln:s","help");


my $help=<<USAGE;
Convert format for alignment

Format aln:
   bl2seq      Bl2seq Blast output
   clustalw    clustalw (.aln) format
   emboss      EMBOSS water and needle format
   fasta       FASTA format
   maf         Multiple Alignment Format
   mase        mase (seaview) format
   mega        MEGA format
   meme        MEME format
   msf         msf (GCG) format
   nexus       Swofford et al NEXUS format
   pfam        Pfam sequence alignment format
   phylip      Felsenstein PHYLIP format
   prodom      prodom (protein domain) format
   psi         PSI-BLAST format
   selex       selex (hmmer) format
   stockholm   stockholm format


Run: perl AlignIO.pl -informat clustalw -outformat fasta -inaln fbox.aln -outaln fbox.fasta

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $in =Bio::AlignIO->new(
                      -file   => $opt{inaln},
                      -format => $opt{informat}
                     );
my $out=Bio::AlignIO->new(
                      -file   => ">$opt{outaln}",
                      -format => $opt{outformat}
                     );

while(my $aln =$in->next_aln()){
     $out->write_aln($aln);
}



