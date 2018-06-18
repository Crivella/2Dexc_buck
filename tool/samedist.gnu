#!/usr/bin/gnuplot

set term pdfcairo enhanced
set encoding utf8
set termoption dashed

set output "plot.pdf"
set fit quiet

FILE="AlN_SAVE.dat"

stats FILE u 2:4 nooutput

xmin = STATS_min_x
xmax = STATS_max_x

f(x) = a + b * x + c * x**2
g(x) = a1 + b1 * x

r=0.07

fit f(x) FILE u 2:4 via a,b,c
set key left top
set title sprintf("f(x)=(%g) + (%g)*x + (%g)*x^2",a,b,c) textcolor rgb "red"

set xrange [xmin:xmax]
set yrange [*:*]

set xlabel "Buckling (bohr)"
set ylabel "[Total energy - a](Ry)"

plot \
	FILE u 2:($4-a) w p pt 7 ps 0.3 t "E_{tot}", \
	f(x)-a w l t "parabolic fit"



fit f(x) FILE u 2:9 via a,b,c
fit [0:r] g(x) FILE u 2:9 via a1,b1
set key right top
set title sprintf("f(x)=(%g) + (%g)*x + (%g)*x^2",a,b,c) textcolor rgb "red"

set xrange [xmin:xmax]
set yrange [*:*]

set xlabel "Buckling (bohr)"
set ylabel "[Γ_{cb} - a](eV)"

set label 1 at graph 0.03, graph 0.2 sprintf("g(x) = (%g) + (%g)*x", a1, b1) textcolor rgb "blue"

plot \
	FILE u 2:($9-a) w p pt 7 ps 0.3 t "Γ_{cb}", \
	f(x)-a w l t "parabolic fit", \
	[0:r] g(x)-a w l t "linear fit"




fit f(x) FILE u 2:8 via a,b,c
fit [0:r] g(x) FILE u 2:8 via a1,b1
set key left top
set title sprintf("f(x)=(%g) + (%g)*x + (%g)*x^2",a,b,c) textcolor rgb "red"

set xrange [xmin:xmax]
set yrange [*:*]

set xlabel "Buckling (bohr)"
set ylabel "[Γ_{vb} - a](eV)"

set label 1 at graph 0.03, graph 0.7 sprintf("g(x) = (%g) + (%g)*x", a1, b1) textcolor rgb "blue"

plot \
	FILE u 2:($8-a) w p pt 7 ps 0.3 t "Γ_{vb}", \
	f(x)-a w l t "parabolic fit", \
	[0:r] g(x)-a w l t "linear fit"






