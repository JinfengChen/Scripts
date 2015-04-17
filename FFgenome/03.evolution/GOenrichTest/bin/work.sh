perl GOlevel.pl -obo ../input/gene_ontology.obo -blast2go ../input/list.txt -level 3 > log &
 perl GOlevel.pl -obo ../input/gene_ontology.obo -blast2go ../input/OB.tandem.GOenrichmentTest.txt -level 3 
> log &

perl Chisq2xXTestR.pl -table FF2otherPfamily4test.txt -number 4 > FF2otherPfamily4test.log 2> log2 &


perl Chisq2xXTestR.pl --table noncollinearTEST --number 2 > noncollinearTEST.log 2> log2 &

perl Chisq2xXTestR.pl --table manaulcheck.noncollinearTEST --number 2 > manaulcheck.noncollinearTEST.test
perl Chisq2xXTestR.pl --table manaulcheck.code1.noncollinearTEST --number 2 > manaulcheck.code1.noncollinearTEST.test



