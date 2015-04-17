#set terminal svg size 600,400 dynamic enhanced fname 'arial'  fsize 10 mousing name "boxplot_1" butt solid 
set term pdfcairo font "Times_New_Roman,8"
set output 'boxplot.pdf'
#set border 2 front linetype -1 linewidth 1.000
set border 2
set boxwidth 0.5 absolute
set style fill   solid 0.25 border lt -1
unset key
set pointsize 0.5
set style data boxplot
set xtics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set xtics  norangelimit
set xtics   ("A" 1.00000, "B" 2.00000)
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set yrange [ 0.00000 : 100.000 ] noreverse nowriteback
plot 'silver.dat' using (1):2, '' using (2):(5*$3)
