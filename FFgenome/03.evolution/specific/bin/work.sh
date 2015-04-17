perl specific.pl -query OB -target OS,SB,BR > log &
perl specific.pl -query OS -target OB,SB,BR > log &
perl specific.pl -query SB -target OB,OS,BR > log &
perl specific.pl -query BR -target OB,SB,OS > log &
