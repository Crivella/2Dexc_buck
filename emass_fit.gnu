#!/usr/bin/gnuplot

set term pdfcairo enhanced
set encoding utf8
set termoption dashed

set output OUTNAME
set print "app.dat"
set fit quiet

stats FILE u 1:5 nooutput

xmin = STATS_min_x
xmax = STATS_max_x

ymin = STATS_min_y
ymax = STATS_max_y
app1 = STATS_pos_min_y
app2 = STATS_pos_max_y

if ( abs( app1 - xmin) < abs( app2 -xmin)) {
	a = ymin
} else {
	a = ymax
}

d=1E-3
f(x) = a + b * x + c * x**2 + d * x**3

fit f(x) FILE u 1:(column(NBND)) via a,b,c,d

print c

set xrange [xmin:xmax]
set yrange [*:*]

set title sprintf("f(x) = %.2f + %.2f *x + %.2f * x^2 + %.2f *x^3",a,b,c,d)
plot \
	FILE u 1:(column(NBND)) w p pt 7 ps 0.3 title "DFT data", \
	f(x) w l title "cubic fit data"





