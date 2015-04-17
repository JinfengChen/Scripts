echo "classify gene into different type of homolog"
perl TypeGene.pl --dir ../input/TIGR6_type > TIGR6.gene.type
perl TypeGene.pl --dir ../input/OBa_type > OBa.gene.type
 
echo "sum gene movement"
perl SumType.pl --list genemovement.txt --type TIGR6.gene.type > genemovement.summary

perl SumType.pl --list ../input/sum/tigr.all.pep.final.pep_pkcat --type TIGR6.gene.type --col 2 > TIGR6.PK.summary
perl SumType.pl --list ../input/sum/OBa.final.pep_pkcat --type OBa.gene.type --col 2 > OBa.PK.summary
perl SumType.pl --list ../input/sum/tigr.all.pep.final.pep_tf_family --type TIGR6.gene.type > TIGR6.TF.summary
perl SumType.pl --list ../input/sum/OBa.final.pep_tf_family --type OBa.gene.type > OBa.TF.summary
