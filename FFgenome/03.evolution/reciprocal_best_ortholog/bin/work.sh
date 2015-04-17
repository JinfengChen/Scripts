#perl reciprocal_best_ortholog.pl ../input/tigr.all.pep.final.fa ../input/glaberrima.final.pep.fa > log 2> log2 &

echo "pairwise ortholog grass"
#perl reciprocal_best_ortholog.pl ../input/OS.pep ../input/OB.pep > log 2> log2 &

perl reciprocal_best_ortholog.pl ../input/OS.pep ../input/BD.pep > log 2> log2
perl reciprocal_best_ortholog.pl ../input/OS.pep ../input/SB.pep > log 2> log2
perl reciprocal_best_ortholog.pl ../input/OS.pep ../input/SI.pep > log 2> log2

perl reciprocal_best_ortholog.pl ../input/OB.pep ../input/BD.pep > log 2> log2
perl reciprocal_best_ortholog.pl ../input/OB.pep ../input/SB.pep > log 2> log2
perl reciprocal_best_ortholog.pl ../input/OB.pep ../input/SI.pep > log 2> log2


