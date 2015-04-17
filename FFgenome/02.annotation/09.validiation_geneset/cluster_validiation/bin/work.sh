cat OBRACH.specific | cut -f 1 > OBRACH.specific.id
perl getidseq.pl -l OBRACH.specific.id -f ../input/OBRACH.pep -o OBRACH.specific.fa
perl runblastall_multicpu.pl -p blastp -i OBRACH.specific.fa -d /home/database/nr/nr > log 2> log2 &
