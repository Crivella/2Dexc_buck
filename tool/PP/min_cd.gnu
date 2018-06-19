#!/usr/bin/gnuplot

if( strstrt( GPVAL_TERMINALS, 'pdfcairo') >0) {
        set term pdfcairo enhanced
} else {
        set term postscript enhanced
}
#set term x11
set encoding utf8
set termoption dashed

PNAME="Etot_vs_cdim.dat"
set print PNAME

PREFIX=system("cat system.sh | grep prefix | cut -d \"=\" -f 2") #"BN_SAVE.dat"
NAME=sprintf("%s_SAVE.dat",PREFIX)

ONAME=sprintf("%s_min_cd.pdf",PREFIX)
set output ONAME

system(sprintf("tools/split_cd.sh %s",NAME))

#LIST="4.64 4.66 4.67 4.674 4.676 4.6765 4.678 4.6785 4.679 4.6795 4.6798 4.6799 4.68"
LIST=system("cat system.sh | grep ALAT_LIST | cut -d \"=\" -f 2 | tr -d '\"' | cut -d \"#\" -f 1")
CLIST="black red blue green magenta cyan orange black red blue green magenta cyan orange black red blue green magenta cyan"
#PLIST="3 3 3 3 3 3 3 7 7 7 7 7 7 7 0 0 0 0 0"
#DLIST="0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 5 5 5"

set xlabel "Buckling (Bohr)"
set ylabel "E_{tot} (Ry)"

p for[i=1:words(LIST)] sprintf("CD_%.5f_%s",word(LIST,i)+0,NAME) u 2:4 w lp lw 2 ps 0.7 \
	lc rgb word(CLIST,i) ti sprintf("celldim %.5f",word(LIST,i)+0), \
  for[i=1:words(LIST)] sprintf("CD_%.5f_%s",word(LIST,i)+0,NAME) u (-$2):4 w lp lw 2 ps 0.7 \
	lc rgb word(CLIST,i) noti 

f(x) = a + b*x + c*x**2
g(x) = e + f*x + g*x**2 + h*x**3

print "#cdim(Bohr)\tEtot(Ry)"
do for[i=1:words(LIST)] {
	CD=word(LIST,i)
	N=sprintf("CD_%.5f_%s",word(LIST,i)+0,NAME)
	fit f(x) N u 2:4 via a,b,c
	app=-b/(2*c)
	print CD,"\t\t", f(app)
	set title sprintf("celldim %.5f",word(LIST,i)+0)
	set key bottom left
	set label 1 at graph 0.2, graph 0.8 tc "red" sprintf("f(x) = %g + %g*x + %g*x^2", a,b,c) font "Verdana,11"
	p N u 2:4 ps 0.5 pt 7  lc rgb "black" ti "data", \
	  N u (-$2):4 ps 0.5 pt 7 lc rgb "black" noti, \
          f(x) lc rgb "red" ti "parabolic fit"
	unset label 1
}

fit f(x) PNAME u 1:2 via a,b,c

set title "E_{tot} vs Cell dimension (minimized for buckling)"
set xlabel "Cell dimension (Bohr)"
set key top right
set label 1 at graph 0.25, graph 0.73 tc "red" sprintf("f(x) = %g + %g*x + %g*x^2", a,b,c) font "Verdana,11"

p PNAME u 1:2 w p pt 7 ps 0.5 lc rgb "black" ti "E_{tot} vs celldim", \
  f(x) lw 2 lc rgb "red" ti "parabolic fit"#, \
  g(x) dt (20,10) lc rgb "blue" ti "cubic"

unset label 1

system("rm CD_*")










