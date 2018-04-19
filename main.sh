#!/bin/bash

source ./ENVIRONMENT_VARIABLES
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND


source ./init.sh
source ./constants.sh

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

	print_str "Running scf calculation buckling=${buck_bohr} bohr" "title" $YELLOW

	IN=${PREFIX}_script.scf_b${buck}.in
	OUT=${PREFIX}_script.scf_b${buck}.out

	#Set k_point for input
	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="12 12 1 1 1 1"

	#Print input for scf
	print_in_pw scf
	#Run command for QE
#	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

	#Extract total energy from output and print it in a two-coloumn file with the buckling
	ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
	echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

	###################################################################################
	#Run band calculation
	echo -e "\n"
	print_str "Running band structure calculation for high-symmetry path of ibrav=$ibrav" "title" $YELLOW
	IN=${PREFIX}_script.nscf_b${buck}.in
	OUT=${PREFIX}_script.nscf_b${buck}.out

	KPT_MODE=""
	KPT_LIST="`cat High_symm/$ibrav.kpt`"

	print_in_pw nscf
#	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

	#Make the band plot 
	BAND_OUT="bands_b${buck}_plotted.dat"
	do_command "qepp_plotband.x $OUT $BAND_OUT" "null"  $BRIGHT_GREEN
	do_command "gnuplot -e FILE='${BAND_OUT}' -e NBND=$nbnd bands.gnu" ""  $BRIGHT_GREEN


	###################################################################################
	#Calcolate effective mass
	#Determine the position of the smallest direct gap and valence and conduction band number
	echo -e "\n"
	print_str "Starting the calcuations for the effective mass" $YELLOW
	IN=""
	APP="$OUT"
	OUT="tmp_gap.out"
	do_command "smallest_gap.x $APP" "io null" $BRIGHT_GREEN

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
	
	MASS_VB=0
	MASS_CB=0
	MASS_VB_1=0
	MASS_CB_1=0
	for app in $AFTER $BEFORE; do
		END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $app | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
		echo -e "${TAB}Running direction $START -> $END"

		IN=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.in
		OUT=${PREFIX}_script.nscf_kpt${N_KPT_MINGAP}to${app}_b${buck}.out

		KPT_MODE="K_POINTS {crystal_b}"
		KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

		print_in_pw nscf
#		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

		BAND_OUT="kpt${N_KPT_MINGAP}to${app}_band.dat"
#		do_command "qepp_plotband.x $OUT $BAND_OUT" "null" $BRIGHT_GREEN
		do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($VB+1)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
		MASS_VB_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
		do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+1)) -e OUTNAME='${BAND_OUT:0: -4}_cb.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
		MASS_CB_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
		rm app.dat

		MIN_LINE=`cat $BAND_OUT | head -n $N_KPT_MINGAP | tail -n 1 |  cut -d$'\t' -f 2`
		E_VB=`echo $MIN_LINE | cut -d" " -f $VB`
		E_VB_1=`echo $MIN_LINE | cut -d" " -f $(($VB-1))`
		DIFF=`echo "$E_VB - $E_VB_1" | bc -l`

		DEGEN=0
		print_str "Looking for degeneracy of the valence band" "sub"
		print_str "  e(vb-1) = $E_VB_1,    e(vb) = $E_VB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			print_str "  Degeneracy found within the 1meV limit.."
			print_str "  running fit of the vb-1 band"
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($VB+0)) -e OUTNAME='${BAND_OUT:0: -4}_vb-1.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
			MASS_VB_1_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
		else
			print_str "  No degeneracy for the valence band"
			print_str "  vb and vb-1 are split more than 1meV"
		fi

		if [[ $DEGEN == "1" ]]; then
			MASS_VB=`echo "$MASS_VB + sqrt(($MASS_VB_APP + $MASS_VB_1_APP)^2)/4" | bc -l`
		else
			MASS_VB=`echo "$MASS_VB + sqrt(($MASS_VB_APP)^2)/2" | bc -l`
		fi

		E_CB=`echo $MIN_LINE | cut -d" " -f $CB`
		E_CB_1=`echo $MIN_LINE | cut -d" " -f $(($CB+1))`
		DIFF=`echo "$E_CB_1 - $E_CB" | bc -l`

		DEGEN=0
		print_str "Looking for degeneracy of the conduction band" "sub"
		print_str "  e(vb-1) = $E_VB_1,    e(vb) = $E_VB  ->  $DIFF"
		if (( `echo "$DIFF<0.001" | bc -l` )); then
			DEGEN=1
			print_str "  Degeneracy found within the 1meV limit.."
			print_str "  running fit of the cb+1 band"
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+2)) -e OUTNAME='${BAND_OUT:0: -4}_cb-1.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
			MASS_CB_1_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
		else
			print_str "  No degeneracy for the conduction band"
			print_str "  vb and vb-1 are split more than 1meV"
		fi

		if [[ $DEGEN == "1" ]]; then
			MASS_CB=`echo "$MASS_CB + sqrt(($MASS_CB_APP + $MASS_CB_1_APP)^2)/4" | bc -l`
		else
			MASS_CB=`echo "$MASS_CB + sqrt(($MASS_CB_APP)^2)/2" | bc -l`
		fi

		#echo "m_vb = $MASS_VB_APP  m_cb = $MASS_CB_APP  m_vb_1 = $MASS_VB_1_APP  m_cb_1 = $MASS_CB_1_APP"
	done

	MU=`echo "(1/$MASS_VB + 1/$MASS_CB)^(-1)" | bc -l`

	#echo "m_h = $MASS_VB     m_e = $MASS_CB   mu = $MU"

	#Performe calculation with dense k-point near minimum before
	END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $BEFORE | tail -n 1 | tr -s " " | cut -d" " -f 2-4`

	####################################################################################
	#Run pw2gw for optical properties
	#Run nscf calculation
	echo -e "\n"
	print_str "Starting the calcuations for the alfa0" "title" $YELLOW

	IN=${PREFIX}_script.nscf-opt_b${buck}.in
	OUT=${PREFIX}_script.nscf-opt_b${buck}.out

	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="30 30 1 0 0 0"

	print_in_pw nscf
#	do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN


	#Run pw2gw
	IN=${PREFIX}_script.pw2gw.in
	OUT=${PREFIX}_script.pw2gw_b${buck}.out

	print_in_pw2gw
#	do_command "$RUN_COMMAND $BIN_DIR/pw2gw_new.x" "date io" $BRIGHT_GREEN

	do_command "eps_average.x epsX.dat epsY.dat" "null" $BRIGHT_GREEN
#	do_command "apply_kk_im.x averaged.dat real_xy.dat 1" "null" $BRIGHT_GREEN

	EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
	VACUUM=`echo "$cdim3 * $alat" | bc -l`
	ALFA0=`echo "$EPS0 * $VACUUM / (4 * $PI)" | bc -l`
	ALFA0=`printf %.5f $ALFA0`

	echo -e "\n${TAB}ALFA0 = $ALFA0"


	####################################################################################
	# Report and final data manipulation
	print_str "Final report for buckling = $buck_bohr (bohr)" "title" $YELLOW
	print_str "m_h = $MASS_VB    m_e = $MASS_CB   red_mass = $MU" "sub" $CYAN
	print_str "alfa_0 = $ALFA0" "sub" $CYAN

	Aex=`echo "1 / $MU" | bc -l`
	Rex=`echo "$MU / 2" | bc -l`
	#interpolate from 'Model_2Dexc.plot
	X=`echo "4*$PI*$ALFA0/$Aex" | bc -l`

#echo $X
	X_N=`echo "($X - 0.101) / 0.001 + 1" | bc -l`
	X_N=`printf %.0f $X_N`
#echo $X_N
	Eb_Rex=`cat model_2Dexc.plot | head -n $X_N | tail -n 1| cut -d$'\t' -f2`
	rex_Aex=`cat model_2Dexc.plot | head -n $X_N | tail -n 1| cut -d$'\t' -f3`

#echo "$Eb_Rex * $Rex * $HA_to_EV"
	Eb=`echo "$Eb_Rex * $Rex * $HA_to_EV" | bc -l`
	rex=`echo "$rex_Aex * $Aex" | bc -l`
	print_str "EXC:  Binding = $Eb   radius = $rex" "sub" $CYAN
done


echo -e "\nEnd at: $(date)"
echo -e "Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
echo -e "\n\n"




exit









