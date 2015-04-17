/opt/training/bin/genscan ../bin/genscan_para/Maize.smat ../input/rice.frags1.fa > rice.frags1.fa.genscan 2> rice.frags1.fa.genscan.log

/opt/training/bin/bgf ../bin/bgf_para/rice.smat ../input/rice.frags1.fa > rice.frags1.fa.bgf 2> rice.frags1.fa.bgf.log

 perl ../bin/predict_convert.pl --predict bgf --perfect --checkcds --final rice --log --verbose ./rice.frags1.fa.bgf ../input/rice.frags1.fa

perl ../bin/predict_convert.pl --predict genscan --perfect --checkcds --final riceg --log --verbose ./rice.frags1.fa.genscan ../input/rice.frags1.fa


