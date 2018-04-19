#!/bin/bash

source ./ENVIRONMENT_VARIABLES
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND
echo "Started at: " `date`

function print_in_pw()
{
cat > $IN << EOF
&control
    calculation       = '$1'
    restart_mode      = 'from_scratch'
    prefix            = '$PREFIX'
    tstress           = .true.
    tprnfor           = .true.
    pseudo_dir        = '$PSEUDO_DIR'
    outdir            = '$TMP_DIR'
    verbosity         = 'high'
    disk_io           = 'minimal'
    wf_collect        = .true.
/
&system
    ibrav             = $ibrav
    celldm(1)         = $alat
    celldm(3)         = $cdim3
    nat               = 2
    ntyp              = 2
    ecutwfc           = $cutoff
    nbnd              = $nbnd
   ! smearing          = 'fermi-dirac'
   ! degauss           = 0.0036749326
/
&electrons
    !diagonalization   = 'cg'
    conv_thr          = 1.0d-8
    mixing_mode       = 'plain'
    mixing_beta       = 0.7d0
/
&ions
    ion_dynamics='bfgs',
    upscale=10
/

ATOMIC_SPECIES
Al  26.982  Al.pz-vbc.UPF
N   14.067  N.pz-vbc.UPF

ATOMIC_POSITIONS {alat}
Al      0.000000000   0.577350300   0.00
N       0.500000000   0.288675100   $buck

$KPT_MODE
$KPT_LIST
EOF
}

function print_in_pw2gw()
{
cat > $IN << EOF
&inputpp
  prefix='$PREFIX',
  outdir='tmp/',
  what='gw',
  qplda=.false.,
  vxcdiag=.false.,
  vkb=.false.,
  Emin=0.0,
  Emax=30.0,
  DeltaE=0.005,
/
EOF
}

export -f print_in_pw
export -f print_in_pw2gw

PREFIX="AlN"

SAVE=Etot_vs_buck.dat
echo -e "# buckling(bohr) Etot(Ry)" > $SAVE

cutoff=30 #100
alat=5.72
cdim3=6
nbnd=20
ibrav=4
#Cycle over different bucklings
for buck in 0.00; do #0.06 0.08 0.10 0.12 0.14
	buck_bohr=`echo $buck*$alat | bc -l`
	buck_bohr=`printf %.4f $buck_bohr`
	echo -e "\t*************************************************************************"
	echo -e "\tRunning scf calculation buckling=${buck_bohr} bohr                       "
	echo -e "\t*************************************************************************"

	IN=${PREFIX}_script.scf_b${buck}.in
	OUT=${PREFIX}_script.scf_b${buck}.out

	#Set k_point for input
	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="12 12 1 1 1 1"

	#Print input for scf
	print_in_pw scf

	#Build commmand for pw.x and execute
	echo -e "\tStart: " `date`
	COMMAND="$RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t  $COMMAND < $IN > $OUT"
	$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#Extract total energy from output and print it in a two-coloumn file with the buckling
	ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
	echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

	###################################################################################
	#Run band calculation
	echo -e "\n"
	echo -e "\t*************************************************************************"
	echo -e "\tRunning band structure calculation for high-symmetry path of ibrav=$ibrav"
	echo -e "\t*************************************************************************"
	IN=${PREFIX}_script.nscf_b${buck}.in
	OUT=${PREFIX}_script.nscf_b${buck}.out

	KPT_MODE=""
	KPT_LIST="`cat High_symm/$ibrav.kpt`"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="$RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t  $COMMAND < $IN > $OUT"
	$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#Make the band plot 
	BAND_OUT="bands_b${buck}_plotted.dat"
	echo -e "\tqepp_plotband.x $OUT $BAND_OUT"
	qepp_plotband.x $OUT $BAND_OUT
	echo -e "\tgnuplot -e \"FILE='$BAND_OUT'\" -e \"NBND=$nbnd\" plot.gnu"
	gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$nbnd" plot.gnu

	###################################################################################
	#Calcolate effective mass
	#Determine the position of the smallest direct gap and valence and conduction band number
	echo -e "\n"
	echo -e "\t*************************************************************************"
	echo -e "\tStarting the calcuations for the effective mass                          "
	echo -e "\t*************************************************************************"
	echo -e "\tsmallest_gap.x $OUT > tmp_gap.out 2>&1"
	smallest_gap.x $OUT > tmp_gap.out 2>&1
	N_KPT_MINGAP=`cat tmp_gap.out | grep "Min gap energy:" -A 1 | tail -n 1 | cut -d# -f2  | tr -dc '0-9,-.'`
	VB=`cat tmp_gap.out | grep "vb = " | cut -d" " -f3 | tr -dc '0-9'`
	CB=`cat tmp_gap.out | grep "vb = " | cut -d" " -f6 | tr -dc '0-9'`
	rm tmp_gap.out

	NUM_KPT=`cat High_symm/$ibrav.kpt | grep -v "#" | wc -l`
	let NUM_KPT-- #Do not count 1st line (K_POINTS {...})
	let NUM_KPT-- #Do not count 2nd line (Number of k-points)
	echo ""
	echo -e "\tMinimal direct gap at kpt number $N_KPT_MINGAP (of $NUM_KPT)"
	echo -e "\tvb = $VB,    cb= $CB"
	let NUM_KPT-- #Do not count last line (By convention the k-path ends with the same point it start with)

	#extract kpt before and after the minimum
	AFTER=$N_KPT_MINGAP
	let AFTER++
	BEFORE=$N_KPT_MINGAP
	let BEFORE--

	if [[ $BEFORE == "0" ]]; then
		BEFORE=$NUM_KPT
	fi

	echo -e "\tMaking nscf calculaiton with dense line from kpt $N_KPT_MINGAP -> $AFTER"
	echo -e "\t                                                 $N_KPT_MINGAP -> $BEFORE"

	#Read starting point from the High_symm k-path
	START=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $N_KPT_MINGAP | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
	
	for app in $AFTER $BEFORE; do
		END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $app | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
		echo -e "\tRunning direction $START -> $END"

		IN=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.in
		OUT=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.out

		KPT_MODE="K_POINTS {crystal_b}"
		KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

		print_in_pw nscf

		echo -e "\tStart: " `date`
		COMMAND="$RUN_COMMAND $BIN_DIR/pw.x"
		echo -e "\t  $COMMAND < $IN > $OUT"
		$COMMAND < $IN > $OUT
		echo -e "\tEnd: " `date`

		BAND_OUT="kpt${N_KPT_MINGAP}to${app}_band.dat"
		echo -e "\tqepp_plotband.x $OUT $BAND_OUT"
		qepp_plotband.x $OUT $BAND_OUT

		#echo -e "\tgnuplot -e \"FILE='$BAND_OUT'\" -e \"NBND=$(($VB+1))\" -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
		echo -e "\tgnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
		echo -e "\tgnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
		gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$(($VB+1))" -e "OUTNAME='${BAND_OUT:0: -4}_vb.pdf'" emass_fit.gnu
		gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$(($CB+1))" -e "OUTNAME='${BAND_OUT:0: -4}_cb.pdf'" emass_fit.gnu

		MIN_LINE=`cat $BAND_OUT | head -n $N_KPT_MINGAP | tail -n 1 |  cut -d$'\t' -f 2`
		E_VB=`echo $MIN_LINE | cut -d" " -f $VB`
		E_VB_1=`echo $MIN_LINE | cut -d" " -f $(($VB-1))`
		DIFF=`echo "$E_VB - $E_VB_1" | bc -l`

		echo -e "\n\t*** Looking for degeneracy of the valence band ***"
		echo -e "\t  e(vb-1) = $E_VB_1,    e(vb) = $E_VB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			echo -e "\t  Degeneracy found within the 1meV limit..."
			echo -e "\t  running fit of the vb-1 band"
			echo -e "\tgnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
			gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$(($VB+0))" -e "OUTNAME='${BAND_OUT:0: -4}_vb-1.pdf'" emass_fit.gnu
		else
			echo -e "\t  No degeneracy for the valence band"
			echo -e "\t  vb and vb-1 are split more than 1meV"
		fi

		E_CB=`echo $MIN_LINE | cut -d" " -f $CB`
		E_CB_1=`echo $MIN_LINE | cut -d" " -f $(($CB+1))`
		DIFF=`echo "$E_CB_1 - $E_CB" | bc -l`
		echo -e "\n\t*** Looking for degeneracy of the conduction band ***"
		echo -e "\t  e(cb) = $E_CB,    e(cb+1) = $E_CB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			echo -e "\t  Degeneracy found within the 1meV limit..."
			echo -e "\t  running fit of the vb-1 band"
			echo -e "\tgnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
			gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$(($CB+2))" -e "OUTNAME='${BAND_OUT:0: -4}_cb-1.pdf'" emass_fit.gnu
		else
			echo -e "\t  No degeneracy for the conduction band"
			echo -e "\t  vb and vb-1 are split more than 1meV"
		fi

	done
	#Performe calculation with dense k-point near minimum before
	END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $BEFORE | tail -n 1 | tr -s " " | cut -d" " -f 2-4`

: '
	IN=${PREFIX}_script.before_b${buck}.in
	OUT=${PREFIX}_script.before_b${buck}.out

	KPT_MODE="K_POINTS {crystal_b}"
	KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="$RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t  $COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#qepp_plotband.x $OUT before_band.dat
' > /dev/null


	####################################################################################
	#Run pw2gw for optical properties
	#Run nscf calculation
	echo -e "\n"
	echo -e "\t*************************************************************************"
	echo -e "\tStarting the calcuations for the alfa0                                   "
	echo -e "\t*************************************************************************"
	IN=${PREFIX}_script.nscf-opt_b${buck}.in
	OUT=${PREFIX}_script.nscf-opt_b${buck}.out

	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="30 30 1 0 0 0"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="$RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t  $COMMAND < $IN > $OUT"
	$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`


	#Run pw2gw
	IN=${PREFIX}_script.pw2gw.in
	OUT=${PREFIX}_script.pw2gw_b${buck}.out

	print_in_pw2gw

	echo -e "\tStart: " `date`
	COMMAND="$RUN_COMMAND $BIN_DIR/pw2gw_new.x"
	echo -e "\t $COMMAND < $IN > $OUT"
	$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	echo -e "\teps_average.x epsX.dat epsY.dat"
	eps_average.x epsX.dat epsY.dat
	echo -e "\tapply_kk_im.x averaged.dat real_xy.dat 1"
	apply_kk_im.x averaged.dat real_xy.dat 1

	EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
	VACUUM=`echo "$cdim3 * $alat" | bc -l`
	ALFA0=`echo "$EPS0 * $VACUUM / (4 * 3.141592)" | bc -l`
	ALFA0=`printf %.5f $ALFA0`

	echo -e"\n\tALFA0 = $ALFA0"
done


echo -e "\nEnd at: $(date)"
echo -e "\n\n"




exit









