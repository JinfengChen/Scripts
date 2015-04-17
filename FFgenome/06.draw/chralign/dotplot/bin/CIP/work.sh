ln -s /share/raid12/chenjinfeng/FFgenome/ortholog/output/blasttable/OB_OS.blasttable ./
ln -s /share/raid12/chenjinfeng/FFgenome/ortholog/output/blasttable/OB_OB.blasttable ./
ln -s /share/raid12/chenjinfeng/FFgenome/ortholog/output/blasttable/OS_OS.blasttable ./
perl ../table2CIPCALP.pl OB_OB.blasttable > ob_ob.blast
perl ../table2CIPCALP.pl OS_OS.blasttable > os_os.blast
perl ../table2CIPCALP.pl OB_OS.blasttable > ob_os.blast
cat ob_ob.blast ob_os.blast os_os.blast > ob
rm *_*
mv ob ob_os.blast

