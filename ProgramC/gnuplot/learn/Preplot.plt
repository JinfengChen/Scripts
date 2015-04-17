set xlabel "Month"
set ylabel "Mass (mm)
set term pdfcairo lw 2 font "Times_New_Roman,8"
set output "precipitation.pdf"
plot "precipitation.dat" u 1:2 w lp pt 5 title "Beijing","precipitation.dat" u 1:3 w lp pt 7 title "Shanghai"
