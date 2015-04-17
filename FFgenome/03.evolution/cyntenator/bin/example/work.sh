../cyntenator -t "((HSX.txt MMX.txt) CFX.txt)" -h phylo HSCFMM.blast "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5)" > test2.txt
../cyntenator -t "((HSX.txt MMX.txt) CFX.txt)" -h blast HSCFMM.blast  > test2.txt 2> log2
