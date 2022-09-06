reset
#set format y "10^{%L}
#set format x "10^{%L}
set terminal postscript  enhanced color "Helvetica,22"
set xlabel "x position (cm)" font "Helvetica,22"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
#
set xrange[0:5]
#set xtics 0,.5,20
set yrange [0:1000]
#
#set logscale y
#set logscale x

set grid
set style line 1 lt 1 lw 4 ps 1 pi 30
set style line 2 lt 3 lw 4 ps 1 pi 30
set style line 3 lt 4 lw 4 ps 1 pi 30
set style line 4 lt 5 lw 4 ps 1 pi 30
set style line 5 lt 6 lw 4 ps 1 pi 30
set style line 6 lt 7 lw 4 ps 1 pi 30

set tmargin 0.1
set rmargin 0.1
set bmargin 0.5
set lmargin 5

set key left top
set title "Rad Hydo Su & Olson Wave"
set ylabel "Temperature (eV)" font "Helvetica,22"
set output "rad_hydro_suolson.eps"
plot "rad_hydro_suolson.dat" u 1:2  w lp ls 3 title "Radiation Temperature", "" u 1:3 w lp ls 2 title "Material Temperature"