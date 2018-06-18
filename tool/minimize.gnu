#!/usr/bin/gnuplot

if( strstrt( GPVAL_TERMINALS, 'pdfcairo') >0) {
        set term pdfcairo enhanced
} else {
        set term postscript enhanced
}
#set term x11
set encoding utf8
set termoption dashed

set output "minimize.pdf"
set print "min.dat"

NAME="BN_SAVE.dat"

LIST="4.64 4.66 4.67 4.674 4.676 4.6765 4.678 4.6785 4.679 4.6795 4.6798 4.6799 4.68"
CLIST="black red blue green magenta cyan orange black red blue green magenta cyan orange black red blue green magenta cyan"
PLIST="3 3 3 3 3 3 3 7 7 7 7 7 7 7 0 0 0 0 0"
DLIST="0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5 5"

set xlabel "Buckling (Bohr)"
set ylabel "E_{tot} (Ry)"

#p for[i=1:words(LIST)] sprintf("%.5f_%s",word(LIST,i)+0,NAME) u 2:4 w lp lw 2 ps 2 \
	pt word(PLIST,i)+0  lc rgb word(CLIST,i) ti sprintf("celldim %.5f",word(LIST,i)+0), \
  for[i=1:words(LIST)] sprintf("%.5f_%s",word(LIST,i)+0,NAME) u (-$2):4 w lp lw 2 ps 2 \
	pt word(PLIST,i)+0 lc rgb word(CLIST,i) noti 

f(x) = a + b*x + c*x**2 + d*x**3

do for[i=1:words(LIST)] {
	CD=word(LIST,i)
	N=sprintf("%.5f_%s",word(LIST,i)+0,NAME)
	#P=word(PLIST,i)
	#C=word(CLIST,i)
	fit f(x) N u 2:4 via a,b,c,d
	app=-2*b/c
	print CD," ", f(app)
	set label 1 at graph 0.2, graph 0.8 sprintf("f(x) = %g + %g*x + %g*x^2 + %g*x^3", a,b,c,d) font "Verdana,9"
	p N u 2:4 ps 0.5 pt 7  lc rgb "black" ti sprintf("celldim %.5f",word(LIST,i)+0), \
	  N u (-$2):4 ps 0.5 pt 7 lc rgb "black" noti, \
          f(x) lc rgb "red"
	unset label 1
}

set title "E_{tot} vs Cell dimension (minimized for buckling)"
fit f(x) "min.dat" u 1:2 via a,b,c,d
set xlabel "Cell dimension (Bohr)"
set label 1 at graph 0.2, graph 0.8 sprintf("f(x) = %g + %g*x + %g*x^2 + %g*x^3", a,b,c,d) font "Verdana,9"
p "min.dat" u 1:2 w p lc rgb "black", \
  f(x) lc rgb "red"

#pt word(PLIST,i)+0 dt(20,word(DLIST,i)+0) lc rgb word(CLIST,i) t

#pause -1



