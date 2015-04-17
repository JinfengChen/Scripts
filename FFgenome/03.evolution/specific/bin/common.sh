formatdb -i OS.specific.fa -p T
blastall -p blastp -i OB.specific.fa -d OS.specific.fa -o OB2OS.blast -e 1e-10 > log &
perl /share/raid12/chenjinfeng/tools/bin/blast_parser.pl -nohead OB2OS.blast > OB2OS.blasttable
cut -f 1 OB2OS.blasttable | sort | uniq | wc -l
cut -f 5 OB2OS.blasttable | sort | uniq | wc -l

OS:286
OB:292
