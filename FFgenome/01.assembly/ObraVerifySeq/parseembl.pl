use Bio::SeqIO;

my $infile=shift;

my $in =Bio::SeqIO->new(-file => 'Hd1ffModify_d.embl',
                        -format => 'EMBL');

my $out = Bio::SeqIO->new (-file => '>hd1ff.fas',
                           -format => 'fasta');


while (my $seq=$in->next_seq()){
      $out->write_seq($seq);
}
