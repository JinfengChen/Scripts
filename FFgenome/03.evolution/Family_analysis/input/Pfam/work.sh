hmmfetch -o NB-ARC.hmm Pfam-A.hmm NB-ARC
hmmfetch -o LRR.hmm -f Pfam-A.hmm LRR.list
hmmfetch -o Pkinase.hmm Pfam-A.hmm Pkinase

/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRR_sd22_2 > LRR_sd22_2.hmm
/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRR_bac_2 >> LRR_sd22_2.hmm
/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRR_RI_2 >> LRR_sd22_2.hmm
/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRRcap_2 >> LRR_sd22_2.hmm
/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRR_typ_2 >> LRR_sd22_2.hmm
/home/biosoftware/hmmer-2.3.2/bin/hmmfetch smart.HMMs LRR_CC_2 >> LRR_sd22_2.hmm
mv LRR_sd22_2.hmm LRR.smart.hmm 
