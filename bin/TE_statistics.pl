#! /usr/local/bin/perl  -w
my $genome_length=shift;
my $gff=shift;
`perl /share/raid006/pengchunfang/GACP/GACP-7.0/03.repeat_finding/denovo_predict/bin/stat_TE2.pl   --gff  $gff --rank type > $gff.TE.subtype`;
`perl /share/raid006/pengchunfang/GACP/GACP-7.0/03.repeat_finding/denovo_predict/bin/after_stat_TE2.pl  $gff.TE.subtype $genome_length`;
`perl /share/raid006/pengchunfang/GACP/GACP-7.0/03.repeat_finding/denovo_predict/bin/statistic_repeat.pl   $gff   $genome_length `;


