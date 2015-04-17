use Bio::SeqIO;

my $infile=shift;

my $in =Bio::SeqIO->new(-file => '20091109genbank.gb',
                        -format => 'genbank');

my $out = Bio::SeqIO->new (-file => '>20091109.fas',
                           -format => 'fasta');


while (my $seq=$in->next_seq()){
      if ($seq->species->binomial =~m/Oryza brachyantha/i){
         $out->write_seq($seq);

      }

}
