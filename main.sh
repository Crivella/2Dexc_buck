#!/bin/bash

source ./ENVIRONMENT_VARIABLES
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND


export ANSII_RESET="\x1B[0m"
export ANSII_F_BLACK="\x1B[30m"
export ANSII_F_RED="\x1B[31m"
export ANSII_F_GREEN="\x1B[32m"
export ANSII_F_YELLOW="\x1B[33m"
export ANSII_F_BLUE="\x1B[34m"
export ANSII_F_MAGENTA="\x1B[35m"
export ANSII_F_CYAN="\x1B[36m"
export ANSII_F_WHITE="\x1B[37m"

export ANSII_F_BRIGHT_RED="\x1B[38;2;255;0;0m"
export ANSII_F_BRIGHT_GREEN="\x1B[38;2;0;255;0m"
export ANSII_F_BRIGHT_BLUE="\x1B[38;2;0;0;255m"
export ANSII_F_BRIGHT_CYAN="\x1B[38;2;0;255;255m"

function set_tab()
{
	TAB=""
	C=0
	while (( $C < $1 )); do
		let C++
		TAB="${TAB}\t"
	done
}

function print_str()
{
	if [[ `echo $2 | grep -i "title"` != "" ]]; then
		echo -e "${TAB}*************************************************************************"
	fi

	if [[ `echo $2 | grep -i "sub"` != "" ]]; then
		echo -e "${TAB}${ANSII_F_RED}***${ANSII_RESET} ${3}${1} ${ANSII_F_RED}***${ANSII_RESET}"
	else
		echo -e "${TAB}${3}${1}${ANSII_RESET}"
	fi

	if [[ `echo $2 | grep -i "title"` != "" ]]; then
		echo -e "${TAB}*************************************************************************"
	fi
}

function do_command()
{
	COMMAND="$1"
	if [[ `echo $2 | grep -i "date"` != "" ]]; then
		echo -e "${TAB}Start: " `date` "   Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
	fi
	
	if [[ `echo $2 | grep -i "io"` != "" ]]; then
		if [[ ${IN} != "" ]]; then
			echo -e "${TAB}  ${3}${COMMAND} < ${IN} > ${OUT}${ANSII_RESET}"
			${COMMAND} < ${IN} > ${OUT}
		else
			echo -e "${TAB}  ${3}${COMMAND} > ${OUT}${ANSII_RESET}"
			${COMMAND} > ${OUT} 2>&1
		fi
	else
		if [[ `echo $2 | grep -i "null"` != "" ]]; then
			echo -e "${TAB}  ${3}${COMMAND} > /dev/null 2>&1${ANSII_RESET}"
			$COMMAND  > /dev/null 2>&1
		else
			echo -e "${TAB}  ${3}${COMMAND}${ANSII_RESET}"
			$COMMAND
		fi
	fi

	if [[ `echo $2 | grep "date"` != "" ]]; then
		echo -e "${TAB}End: " `date` "   Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
	fi
}

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
export -f do_command
export -f print_str

set_tab 0
echo "Started at: " `date`
START_TIME=`date +%s.%N`

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
	set_tab 1
	buck_bohr=`echo $buck*$alat | bc -l`
	buck_bohr=`printf %.4f $buck_bohr`

	print_str "Running scf calculation buckling=${buck_bohr} bohr" "title" $ANSII_F_YELLOW

	IN=${PREFIX}_script.scf_b${buck}.in
	OUT=${PREFIX}_script.scf_b${buck}.out

	#Set k_point for input
	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="12 12 1 1 1 1"

	#Print input for scf
	print_in_pw scf
	#Run command for QE
	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $ANSII_F_BRIGHT_GREEN

	#Extract total energy from output and print it in a two-coloumn file with the buckling
	ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
	echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

	###################################################################################
	#Run band calculation
	echo -e "\n"
	print_str "Running band structure calculation for high-symmetry path of ibrav=$ibrav" "title" $ANSII_F_YELLOW
	IN=${PREFIX}_script.nscf_b${buck}.in
	OUT=${PREFIX}_script.nscf_b${buck}.out

	KPT_MODE=""
	KPT_LIST="`cat High_symm/$ibrav.kpt`"

	print_in_pw nscf
	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $ANSII_F_BRIGHT_GREEN

	#Make the band plot 
	BAND_OUT="bands_b${buck}_plotted.dat"
	do_command "qepp_plotband.x $OUT $BAND_OUT" "null"  $ANSII_F_BRIGHT_GREEN
	do_command "gnuplot -e FILE='${BAND_OUT}' -e NBND=$nbnd plot.gnu" ""  $ANSII_F_BRIGHT_GREEN


	###################################################################################
	#Calcolate effective mass
	#Determine the position of the smallest direct gap and valence and conduction band number
	echo -e "\n"
	print_str "Starting the calcuations for the effective mass" $ANSII_F_YELLOW
	IN=""
	APP="$OUT"
	OUT="tmp_gap.out"
	do_command "smallest_gap.x $APP" "io null" $ANSII_F_BRIGHT_GREEN

	N_KPT_MINGAP=`cat tmp_gap.out | grep "Min gap energy:" -A 1 | tail -n 1 | cut -d# -f2  | tr -dc '0-9,-.'`
	VB=`cat tmp_gap.out | grep "vb = " | cut -d" " -f3 | tr -dc '0-9'`
	CB=`cat tmp_gap.out | grep "vb = " | cut -d" " -f6 | tr -dc '0-9'`
	rm tmp_gap.out

	NUM_KPT=`cat High_symm/$ibrav.kpt | grep -v "#" | wc -l`
	let NUM_KPT-- #Do not count 1st line (K_POINTS {...})
	let NUM_KPT-- #Do not count 2nd line (Number of k-points)
	echo ""
	print_str "Minimal direct gap at kpt number $N_KPT_MINGAP (of $NUM_KPT)"
	print_str "vb = $VB,    cb= $CB"
	let NUM_KPT-- #Do not count last line (By convention the k-path ends with the same point it start with)

	#extract kpt before and after the minimum
	AFTER=$N_KPT_MINGAP
	let AFTER++
	BEFORE=$N_KPT_MINGAP
	let BEFORE--

	if [[ $BEFORE == "0" ]]; then
		BEFORE=$NUM_KPT
	fi

	print_str "Making nscf calculaiton with dense line from kpt $N_KPT_MINGAP -> $AFTER"
	print_str "                                                 $N_KPT_MINGAP -> $BEFORE"

	#Read starting point from the High_symm k-path
	START=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $N_KPT_MINGAP | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
	
	for app in $AFTER $BEFORE; do
		END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $app | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
		echo -e "${TAB}Running direction $START -> $END"

		IN=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.in
		OUT=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.out

		KPT_MODE="K_POINTS {crystal_b}"
		KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

		print_in_pw nscf
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $ANSII_F_BRIGHT_GREEN

		BAND_OUT="kpt${N_KPT_MINGAP}to${app}_band.dat"
		do_command "qepp_plotband.x $OUT $BAND_OUT" "null" $ANSII_F_BRIGHT_GREEN
		do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($VB+1)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $ANSII_F_BRIGHT_GREEN
		do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+1)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $ANSII_F_BRIGHT_GREEN

		MIN_LINE=`cat $BAND_OUT | head -n $N_KPT_MINGAP | tail -n 1 |  cut -d$'\t' -f 2`
		E_VB=`echo $MIN_LINE | cut -d" " -f $VB`
		E_VB_1=`echo $MIN_LINE | cut -d" " -f $(($VB-1))`
		DIFF=`echo "$E_VB - $E_VB_1" | bc -l`

		print_str "Looking for degeneracy of the valence band" "sub"
		print_str "  e(vb-1) = $E_VB_1,    e(vb) = $E_VB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			print_str "  Degeneracy found within the 1meV limit.."
			print_str "  running fit of the vb-1 band"
			print_str "gnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($VB+0)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $ANSII_F_BRIGHT_GREEN
		else
			print_str "  No degeneracy for the valence band"
			print_str "  vb and vb-1 are split more than 1meV"
		fi

		E_CB=`echo $MIN_LINE | cut -d" " -f $CB`
		E_CB_1=`echo $MIN_LINE | cut -d" " -f $(($CB+1))`
		DIFF=`echo "$E_CB_1 - $E_CB" | bc -l`
		print_str "Looking for degeneracy of the conduction band" "sub"
		print_str "  e(vb-1) = $E_VB_1,    e(vb) = $E_VB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			print_str "  Degeneracy found within the 1meV limit.."
			print_str "  running fit of the vb-1 band"
			print_str "gnuplot -e \"OUTNAME='${BAND_OUT:0: -4}_vb.pdf'\" emass_fit.gnu"
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+2)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $ANSII_F_BRIGHT_GREEN
		else
			print_str "  No degeneracy for the conduction band"
			print_str "  vb and vb-1 are split more than 1meV"
		fi

	done

	#Performe calculation with dense k-point near minimum before
	END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $BEFORE | tail -n 1 | tr -s " " | cut -d" " -f 2-4`

	####################################################################################
	#Run pw2gw for optical properties
	#Run nscf calculation
	echo -e "\n"
	print_str "Starting the calcuations for the alfa0" "title" $ANSII_F_YELLOW

	IN=${PREFIX}_script.nscf-opt_b${buck}.in
	OUT=${PREFIX}_script.nscf-opt_b${buck}.out

	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="30 30 1 0 0 0"

	print_in_pw nscf
	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $ANSII_F_BRIGHT_GREEN


	#Run pw2gw
	IN=${PREFIX}_script.pw2gw.in
	OUT=${PREFIX}_script.pw2gw_b${buck}.out

	print_in_pw2gw
	do_command "$RUN_COMMAND $BIN_DIR/pw2gw_new.x" "date io" $ANSII_F_BRIGHT_GREEN

	do_command "eps_average.x epsX.dat epsY.dat" "null" $ANSII_F_BRIGHT_GREEN
	do_command "apply_kk_im.x averaged.dat real_xy.dat 1" "null" $ANSII_F_BRIGHT_GREEN

	EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
	VACUUM=`echo "$cdim3 * $alat" | bc -l`
	ALFA0=`echo "$EPS0 * $VACUUM / (4 * 3.141592)" | bc -l`
	ALFA0=`printf %.5f $ALFA0`

	echo -e "\n${TAB}ALFA0 = $ALFA0"
done


echo -e "\nEnd at: $(date)"
echo -e "Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
echo -e "\n\n"




exit









