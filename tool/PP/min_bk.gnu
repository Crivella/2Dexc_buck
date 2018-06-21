#!/usr/bin/gnuplot

if( strstrt( GPVAL_TERMINALS, 'pdfcairo') >0) {
        set term pdfcairo enhanced
} else {
        set term postscript enhanced
}
#set term x11
set encoding utf8
set termoption dashed

PNAME="Etot_vs_bk.dat"
set print PNAME

PREFIX=system("cat system.sh | grep prefix | cut -d \"=\" -f 2") #"BN_SAVE.dat"
NAME=sprintf("%s_SAVE.dat",PREFIX)

ONAME=sprintf("%s_min_bk.pdf",PREFIX)
set output ONAME

system(sprintf("tool/PP/split_bk.sh %s",NAME))

#LIST="0.0000 0.0100 0.0200 0.0300 0.0500 0.0700 0.0900"
LIST=system("cat system.sh | grep BUCKLING_LIST | cut -d \"=\" -f 2 | tr -d '\"' | cut -d \"#\" -f 1")
CLIST="black red blue green magenta cyan orange black red blue green magenta cyan orange black red blue green magenta cyan"

ALAT_0=PREFIX=system("cat system.sh | grep ALAT_0= | cut -d \"=\" -f 2 | cut -d \"#\" -f 1")

set xlabel "Cell dimension (Bohr)"
set ylabel "E_{tot} (Ry)"

set key outside font "Verdana,7"
p for[i=1:words(LIST)] sprintf("BK_%.4f_%s",word(LIST,i)+0,NAME) u 1:4 w lp lw 2 ps 0.5 \
	lc rgb word(CLIST,i) ti sprintf("buckling %.4f(alat)",word(LIST,i)+0)

set key inside font "Verdana,12"

f(x) = a + b*x + c*x**2
g(x) = e + d*x

print "#buck(alat)\tEtot(Ry)"
do for[i=1:words(LIST)] {
	BUCK=word(LIST,i)
	N=sprintf("BK_%.4f_%s",word(LIST,i)+0,NAME)
	fit f(x) N u 1:4 via a,b,c
	fit g(x) N u 1:3 via e,d
	app=-b/(2*c)
	dist=g(app)
	print BUCK,"\t\t", f(app)
	set label 1 at graph 0.25, graph 0.8 tc "red" sprintf("f(x) = %g + %g*x + %g*x^2", a,b,c) font "Verdana,11"
	set label 2 at graph 0.35, graph 0.7 tc "blue" sprintf("dist = %g (bohr)", dist) font "Verdana,11"
	set label 3 at graph 0.35, graph 0.6 tc "green" sprintf("min = %g", app) font "Verdana,11"
	p N u 1:4 ps 0.5 pt 7  lc rgb "black" ti sprintf("buckling %.5f(bohr)",(word(LIST,i)+0)*ALAT_0), \
          f(x) lc rgb "red" ti "parabolic fit"
	unset label 1
	unset label 2
	unset label 3
}

a=0; b=0; c=0
fit f(x) PNAME u ($1*4.68):2 via a,b,c

set title "E_{tot} vs buckling (minimized for cell dimension)"
set xlabel "Buckling (Bohr)"
set ylabel "E_{tot} - a (Ry)"
set key top left
set label 1 at graph 0.2, graph 0.8 tc "red" sprintf("f(x) = %g + %g*x + %g*x^2", a,b,c) font "Verdana,11"

p PNAME u ($1*4.68):($2-a) w p pt 7 ps 0.5 lc rgb "black" ti "E_{tot} vs buckling", \
  f(x)-a lc rgb "red" ti "parabolic fit"

unset label 1

system("rm BK_*")









