perl modifyfsa.pl -i ./fsa/1p3a.fsa -l DACP -protein Dhydrolipoyllysine-residue acetyltransferase > ./gbf/DACP.fsa
perl modifyfsa.pl -i ./fsa/1p8.fsa -l SIGT -protein Signal transducer/ two-component sensor molecule > ./gbf/SIGT.fsa
perl modifyfsa.pl -i ./fsa/2p1.fsa -l ADTY -protein ATP-dependent transporter YFL028C > ./gbf/ADTY.fsa
perl modifyfsa.pl -i ./fsa/3p3.fsa -l PP2C -protein Catalytic/ protein phosphatase type 2C > ./gbf/PP2C.fsa
perl modifyfsa.pl -i ./fsa/8p1.fsa -l SPS1 -protein Sucrose-phosphate synthase 1 > ./gbf/SPS1.fsa
perl modifyfsa.pl -i ./fsa/10p1.fsa -l UPL -protein Ubiquitin-protein ligase > ./gbf/UPL.fsa
perl modifyfsa.pl -i ./fsa/10p3.fsa -l TIFIIF -protein Transcription initiation factor IIF, alpha subunit > ./gbf/TIFIIF.fsa
perl modifyfsa.pl -i ./fsa/11p1.fsa -l TRAN -protein Transferase, transferring glycosyl groups  > ./gbf/TRAN.fsa
perl modifyfsa.pl -i ./fsa/12p3.fsa -l MDEH -protein Malate dehydrogenase, glyoxysomal precursor > ./gbf/MDEH.fsa
tbl2asn -t foxtail.sbt -p ./gbf -a l -V b
