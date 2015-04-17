# set terminal svg size 600,400 dynamic enhanced fname 'arial'  fsize 10 mousing name "boxplot_4" butt solid 
# set output 'boxplot.4.svg'
set term pdfcairo font "Times_New_Roman,8"
set output "boxplot_2.pdf"
set border 2 front linetype -1 linewidth 1.000
set boxwidth 0.5 absolute
set style fill   solid 0.25 border lt -1
unset key
set pointsize 0.5
set style data boxplot
set xtics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set xtics  norangelimit
set xtics   ("Coal" 1.00000, "Gas" 2.00000, "Hydroelectric" 3.00000, "Nuclear" 4.00000, "Oil" 5.00000, "Renewable" 6.00000)
set ytics border in scale 1,0.5 nomirror norotate  offset character 0, 0, 0 autojustify
set title "Distribution of energy usage of the continents, grouped by type of energy source\n" 
set ylabel "Billion Tons of Oil Equivalent" 
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto sorted
t(x) = x/1.e6
filter(col, factor_col, level) = (strcol(factor_col) eq word(factors, level)) ? t(column(col)) : 1/0
factors = "Coal Gas Hydroelectric Nuclear Oil Renewable"
n_f = 6
i = 6
GPFUN_t = "t(x) = x/1.e6"
GPFUN_filter = "filter(col, factor_col, level) = (strcol(factor_col) eq word(factors, level)) ? t(column(col)) : 1/0"
plot for [i=1:n_f] 'energy_circles.dat' using (i):(filter(8, 4, i))
set output

