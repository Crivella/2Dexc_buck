#!/usr/bin/gnuplot

if( strstrt( GPVAL_TERMINALS, 'pdfcairo') >0) {
	set term pdfcairo enhanced
} else {
	set term postscript enhanced
}
set encoding utf8
set termoption dashed

set output OUTNAME

stats FILE u 1:5 nooutput

xmin = STATS_min_x
xmax = STATS_max_x

set xrange [xmin:xmax]
set yrange [-16:8]

set arrow 1 from xmin,0 to xmax,0 nohead dt(10,5) lc rgb "grey"

plot \
	for[i=2:NBND] FILE u 1:i w l notitle


