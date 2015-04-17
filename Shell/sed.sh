echo "sed also output these line that not match, perl will not"
grep "mRNA" tigr.all.final.gff | sed 's/.*ID=\(.*\)\;Name.*/\1/' > tigr.all.final.id
grep "mRNA" tigr.all.final.gff | perl -e'while (my $line=<>){if($line=~/ID\=(.*?);Name/){print $1,"\n"}}' > tigr.all.final.id
echo "add -n and p to sed can only output matched line"
sed -n 's/ab\(.*\)dd/\1/p' test.txt
test.txt
abcdadd
adba
dfdfd

