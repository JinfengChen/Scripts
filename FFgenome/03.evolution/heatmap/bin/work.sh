perl headmapPipe.pl -chr chr12 -riceTE ../input/IRGSP.build5.RepeatMasker.out.gff.chr -ffTE ../input/OBa.all.fa.RepeatMasker.out.gff.chr -ffGene ../input/OBa.all.gff.chr -riceGene ../input/RAP3.gff3.nr.gff.chr > log 2> log2 &

perl remove_redundance.pl ../input/chr12/FF.RT.gff ../input/chr12/FF.DNA.gff ../input/chr12/FF.Gene.gff > log 2> log2 &


perl distri_data_pre.pl ../input/fflen distri.gff.nr.out > log 2> log2 &

perl get_sub_gff.pl -list ../input/ffpara.gene -gff ../input/chr12/FF.Gene.gff -output ffpara.gene.gff

perl get_paralog_exons.pl ffpara.gene.gff

perl get_introns_exons.pl ../input/chr12/FF.Gene.gff

perl density_data_pre.pl ../input/fflen ../input/chr12/FF.GYPSY.gff ../input/chr12/FF.COPIA.gff ../input/chr12/FF.CACTA.gff ../input/chr12/FF.MITE.gff gene_exons.gff paralog_exons.gff > log 2> log2 &

perl Draw_HeatMap.v2.pl A12.data.distri A12.data.densi > log 2> log2 &

/share/raid12/chenjinfeng/tools/draw/svg2xxx_release/svg2xxx -t png -m 400 A12.HeatMap.svg

run chr11 same as chr12

perl connect.pl -table ../input/ob_os.blast -gff1 ../input/chr11/FF.Gene.gff -gff2 ../input/chr12/FF.Gene.gff > log 2> log2 &

perl Draw_HeatMap.v3.pl A12.data.distri A12.data.densi A11.data.distri A11.data.densi all.gene.position connect.txt > log 2> log2 &

