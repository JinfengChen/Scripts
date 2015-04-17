perl /nas/GAG_02/minjiumeng/GACP-8.0/01.gene_finding/cDNA-EST-map-genome/bin/est_map_genome.pl --dbcut 10 --tophit 3 --cpu 50 --run qsub --identity 0.95 --alignrate 0.9 Os.seq.all.change.fa /ifs1/GAG/annotation/minjiumeng/FF_3.0/00.data/Oryza_brachyantha.genome.super_scaffold.v1.0.fa &
perl /nas/GAG_02/minjiumeng/GACP-8.0/01.gene_finding/cDNA-EST-map-genome/bin/sim4_to_gff3.pl --identity 0.9 --alignrate 0.8 --tophit 1 Os.seq.all.change.fa.blat.sim4 >Os.seq.all.change.fa.blat.sim4.gff
nohup perl /nas/GAG_02/minjiumeng/GACP-8.0/01.gene_finding/cDNA-EST-map-genome/bin/pasa_alignment_assembly_by_cluster.pl /nas/GAG_02/minjiumeng/GACP-8.0/01.gene_finding/cDNA-EST-map-genome/software/PASA/pasa_cpp/pasa Os.seq.all.change.fa.blat.sim4.gff >Os.seq.all.change.fa.blat.sim4.pasa.gff &

