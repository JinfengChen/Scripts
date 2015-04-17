perl dotplotRegion.pl --refseq Moc1.fasta --refgenegff Moc1.gff --qryseq Moc1_Scaffold000018_3823279_4020080.fasta --qrygenegff Moc1_Scaffold000018_3823279_4020080.gff --project test

perl dotplotRegion.pl --refseq OB_chr12_13158000_13450000.fasta --qryseq OS_chr12_24586000_25020000.fasta --refgenegff OB_chr12_13158000_13450000.gff --qrygenegff OS_chr12_24586000_25020000.gff --compare OB_chr12_13158000_13450000VSOS_chr12_24586000_250200004ACT --type ACT --project test > log 2> log2 &

perl drawACT.pl --refseq OB_chr12_13232000_13276000.fasta --qryseq OS_chr12_24746000_24792000.fasta --refgenegff OB_chr12_13232000_13276000.gff --qrygenegff OS_chr12_24746000_24792000.gff --compare OB_chr12_13232000_13276000VSOS_chr12_24746000_247920004ACT --type ACT --project test


