#!/usr/bin/gnuplot

if( strstrt( GPVAL_TERMINALS, 'pdfcairo') >0) {
        set term pdfcairo enhanced
} else {
        set term postscript enhanced
}
#set term x11
set encoding utf8
set termoption dashed

set output "inplane.pdf"

NAME="BN_SAVE.dat"

f(x) = a + b*x + c*x**2
#g(x) = a + b*x + c*x**2 + d*x**3

system(sprintf("grep '0.00000' %s > 1.tab",NAME))

set xlabel "Cell dimension (Bohr)"
set ylabel "E_{tot} (Ry)"

fit f(x) "1.tab" u 1:4 via a,b,c 

set title "In plane deviations"
set label 1 at graph 0.3, graph 0.8 tc "red" sprintf("f(x) = %g + %g*x + %g*x^2", a,b,c) font "Verdana,9.5"

p "1.tab" u 1:4 ps 0.5 pt 7 lc rgb "black" ti "E_{tot} vs celldim", \
  f(x) lc rgb "red" ti "parabolic fit"

system("rm 1.tab")

#pause -1





