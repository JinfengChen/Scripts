perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:len /share/raid12/chenjinfeng/FFgenome/assmbly/superscaf/bin/bigscaffold/OBa.chr.fa > FF.chr.len
perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:len /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/rice_data/IRGSP.build5 > IRGSP5.chr.len
paste IRGSP5.chr.len FF.chr.len > os_ob.chr
