echo "OBvsOS"
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/OS > log 2> log2 
perl runblastall.pl -p blastp -f T -i ../input/OS -d ../input/OB > log 2> log2 
perl runblastall.pl -p blastp -f T -i ../input/OS -d ../input/OS > log 2> log2 
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/OB > log 2> log2 
echo "OBvsSB"
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/SB > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/SB -d ../input/OB > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/SB -d ../input/SB > log 2> log2
echo "OBvsBR"
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/BR > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/BR -d ../input/OB > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/BR -d ../input/BR > log 2> log2
echo "OSvsSB"
perl runblastall.pl -p blastp -f T -i ../input/OS -d ../input/SB > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/SB -d ../input/OS > log 2> log2
echo "OSvsBR"
perl runblastall.pl -p blastp -f T -i ../input/OS -d ../input/BR > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/BR -d ../input/OS > log 2> log2
echo "SBvsBR"
perl runblastall.pl -p blastp -f T -i ../input/SB -d ../input/BR > log 2> log2
perl runblastall.pl -p blastp -f T -i ../input/BR -d ../input/SB > log 2> log2
echo "SBvsZM"
perl runblastall.pl -p blastp -f T -i ../input/ZM -d ../input/ZM > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/SB -d ../input/ZM > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/ZM -d ../input/SB > log 2> log2 &

echo "OSvsZM"
perl runblastall.pl -p blastp -f T -i ../input/ZM -d ../input/OS > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/OS -d ../input/ZM > log 2> log2 &

echo "OBvsZM"
perl runblastall.pl -p blastp -f T -i ../input/ZM -d ../input/OB > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/ZM > log 2> log2 &

echo "BRvsZM"
perl runblastall.pl -p blastp -f T -i ../input/ZM -d ../input/BR > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/BR -d ../input/ZM > log 2> log2 &

echo "TRvsTR" ## tigr6 
perl runblastall.pl -p blastp -f T -i ../input/TR -d ../input/TR > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/TR -d ../input/OB > log 2> log2 &
perl runblastall.pl -p blastp -f T -i ../input/OB -d ../input/TR > log 2> log2 &

