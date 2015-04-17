g++ -o plot plot.cpp

echo "this gnuplot supports pdfcario"
../gnuplot/bin/gnuplot Preplot.plt

echo "this gnuplot support boxplot"
../gnuplot4.6/bin/gnuplot boxplot.plt
../gnuplot4.6/bin/gnuplot boxplot_2.plt
