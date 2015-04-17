perl /share/raid12/chenjinfeng/FFgenome/evolution/GOboxplot/bin/GOparse.pl -obo /share/raid12/chenjinfeng/FFgenome/evolution/GOboxplot/input/gene_ontology.obo -wego /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/Sorghum_bicolor/Sb.final.pep.iprscan.gene.wego -table /share/raid12/chenjinfeng/FFgenome/ortholog/output/orthoPair/KsFigure/KaKs_calculator/SB-ZM.KaKs.summary.txt -column 5 -level 3 -project SB2ZM_Identity > log 2> log2 &

perl /share/raid12/chenjinfeng/FFgenome/evolution/GOboxplot/bin/boxplot.pl -class BP -table SB2ZM_Identity.GO.result -project SB2ZM_Identity > log 2> log2 &


