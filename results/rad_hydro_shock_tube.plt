reset
#set format y "10^{%L}
#set format x "10^{%L}
set terminal postscript  enhanced color "Helvetica,22"
set xlabel "x position (cm)" font "Helvetica,22"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
#
#set xrange[0.0:3.5]
#set xtics 0,150,1500
#set yrange [1:3.5]
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

set tmargin 0.5
set rmargin 0.5
set bmargin 0.5
set lmargin 7


set key left top
set title "Non-Equilibrium Shock Radiation Vs Material Temperature"
set ylabel "Temperature (eV)" font "Helvetica,22"
set output "rad_hydro_shock_tube_temp.eps"
plot "rad_hydro_shock_so.dat" u 1:2  w lp ls 3 title "Radiation Temperature", "" u 1:3 w lp ls 2 title "Material Temperature"

set key left top
set title "Non-Equilibrium Shock Density"
set ylabel "Density (g/cc)" font "Helvetica,22"
set output "rad_hydro_shock_tube_dens.eps"
plot "rad_hydro_shock_so.dat" u 1:4  w lp ls 3 title "Density"

set key left top
set title "Non-Equilibrium Shock Velocity"
set ylabel "Velocity" font "Helvetica,22"
set output "rad_hydro_shock_tube_vel.eps"
plot "rad_hydro_shock_so.dat" u 1:5  w lp ls 3 title "Velocity"

set key left top
set title "Non-Equilibrium Shock Total Energy"
set ylabel "Velocity" font "Helvetica,22"
set output "rad_hydro_shock_tube_ener.eps"
plot "rad_hydro_shock_so.dat" u 1:6  w lp ls 3 title "Total Energy"

set key left top
set title "Non-Equilibrium Shock Total Pressure"
set ylabel "Pressure" font "Helvetica,22"
set output "rad_hydro_shock_tube_pres.eps"
plot "rad_hydro_shock_so.dat" u 1:7  w lp ls 3 title "Total Pres"