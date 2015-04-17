### must be run at where the sequence is (perl yn00.pl ./)

if (@ARGV < 1){
   print "Enter a directory: './'\n";
   exit ();
}

my $dir=$ARGV[0];
print "Translating and Aligning\n";
`perl /share/raid12/chenjinfeng/FFgenome/evolution/kaks_pairwise/bin/script_yn00/dna2aarunclustal.pl $dir > dna2aa.log 2> dna2aa.log2`;
print "Convert aa align to dna align\n";
`perl /share/raid12/chenjinfeng/FFgenome/evolution/kaks_pairwise/bin/script_yn00/aa2dnaalign.pl $dir > aa2dna.log 2> aa2dna.log2`;
print "Fasta to paml formating\n";
`perl /share/raid12/chenjinfeng/FFgenome/evolution/kaks_pairwise/bin/script_yn00/fasta2paml.pl $dir > fasta2paml.log 2> fasta2paml.log2`;
print "Done\n";
