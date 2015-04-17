echo "prepare cds file"
perl cleanhead.pl --fasta ../input/tigr.all.cds.final.fa --clean ../input/tigr.all.cds.final.clean.fa
perl cleanhead.pl --fasta ../input/Bradi_1.0.cds.final.fa --clean ../input/Bradi_1.0.cds.final.clean.fa
perl cleanhead.pl --fasta ../input/Sb.final.cds --clean ../input/Sb.final.clean.fa
echo "prepare pos file"
perl cleanheadgff.pl -gff ../input/Bradi_1.0.final.gff > gene.pos
perl cleanheadgff.pl -gff ../input/tigr.all.final.gff >> gene.pos
perl cleanheadgff.pl -gff ../input/Sb.final.gff >> gene.pos


perl CloseHomolog.pl --blast ../input/Bradi_1.0.cds.final.clean.fa.blasttable > log 2> log2 &
perl CloseHomolog.pl --blast ../input/Bradi_1.0.cds.final.clean.fa.blasttable --position ../input/gene.pos > log 2> log2 &
mv log global_synteny_table

echo "find colinear gene among species"
perl gene_movement_mk_colinear_gene_list.pl global_synteny_table > colinear_gene_list1
echo "find some non-colinear gene that moved within regions"
perl gene_movement_mk_colinear_gene_list2.pl colinear_gene_list1 
echo "edit select line in gene_movement_mk_colinear_gene_list.pl to do non-colinear gene selection"
perl gene_movement_mk_colinear_gene_list.pl global_synteny_table > Brd.0.non_colinear_genes
echo "find donor sites"
perl gene_movement_donor_sites.pl colinear_genes non_colinear_genes > log 2> log2 &

perl gene_movement_donor_sites.pl colinear_gene_list2 Brd.0.non_colinear_genes > log 2> log2 &

echo "get region data"
perl getsubdata.pl --genegff ../input/Bradi_1.0.final.gff --fasta ../input/Brachypodium_distachyon_Bd21.main_genome.scaffolds.fasta --refgene Bradi1g78730.1 --qrygene Bradi1g48110.1

echo "dotplot region data"
perl dotplotRegion.pl --refseq Brd_Bd1_46693171_46745366.fasta --qryseq Brd_Bd1_74600312_74651997.fasta --refgenegff Brd_Bd1_46693171_46745366.gene.gff --qrygenegff Brd_Bd1_74600312_74651997.gene.gff --project test

echo "get region data and draw"
perl getsubdata.pl --genegff ../input/Bradi_1.0.final.gff --fasta ../input/Brachypodium_distachyon_Bd21.main_genome.scaffolds.fasta --refgene Bradi1g78730.1 --qrygene Bradi1g48110.1 --flanking 10000 --draw

