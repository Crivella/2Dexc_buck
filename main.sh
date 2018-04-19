#!/bin/bash

source ./ENVIRONMENT_VARIABLES
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND


source ./init.sh
source ./constants.sh

TAB_C=0
set_tab $TAB_C
echo "Started at: " `date`
START_TIME=`date +%s.%N`

PREFIX="AlN"

SAVE=Etot_vs_buck.dat
echo -e "# buckling(bohr) Etot(Ry)" > $SAVE

cutoff=30 #100
ALAT=5.72
CDIM3=6
nat=2
ntyp=2
nbnd=20
ibrav=4

ATOM_1="Al  26.982  Al.pz-vbc.UPF"
ATOM_2="N   14.067  N.pz-vbc.UPF"
ATOM_3=""
ATOM_4=""
ATOM_5=""
ATOM_6=""

AP="alat"
POS_1="Al      0.000000000   0.577350300   0.00"
POS_2="N       0.500000000   0.288675100   0.00"
POS_3=""
POS_4=""
POS_5=""
POS_6=""
POS_7=""
POS_8=""
POS_9=""
POS_10=""
POS_11=""
POS_12=""

ALAT_LIST="5.72 5.70 5.74"
BUCKLING_LIST="0.00"

#Cycle over different celldimension
A1_x=`echo $POS_1 | tr -s " " | cut -d" " -f 2`
A1_y=`echo $POS_1 | tr -s " " | cut -d" " -f 3`
A2_x=`echo $POS_2 | tr -s " " | cut -d" " -f 2`
A2_y=`echo $POS_2 | tr -s " " | cut -d" " -f 3`
DIST_alat=`echo "sqrt(($A1_x - $A2_x)^2 + ($A1_y - $A2_y)^2)" | bc -l`
DIST_0=`echo "$ALAT * $DIST_alat" | bc -l`

let TAB_C++
for alat in $ALAT_LIST; do
	cdim3=`echo "$ALAT*$CDIM3 / $alat" | bc -l`
	cdim3=`printf %.4f $cdim3`

#echo "cdim3 = $cdim3"
	#Introduce buckling value to calculate that gives atom distance equal to DIST_0
	#Fixed atom distance calculation
	if (( `echo "$alat <= $ALAT" | bc -l` )); then
		CHECK=0
		for b in $BUCKLING_LIST; do
			DIST=`echo "$alat * sqrt(($A1_x - $A2_x)^2 + ($A1_y - $A2_y)^2 + $b^2)" | bc -l`
#echo "$DIST - $DIST_0"
			if (( `echo "sqrt(($DIST - $DIST_0)^2)<0.001" | bc -l` )); then
				CHECK=1
				break
			fi
		done
		blist="$BUCKLING_LIST"
		if [[ $CHECK == "0" ]]; then
			b=`echo "sqrt($DIST_0 ^ 2 - ($DIST_alat * $alat)^2)" | bc -l`
			b=`printf %.5f $b`
			blist="$BUCKLING_LIST $b"
			print_str "Added buckling value $b to loop to make constant distant calculation"
		fi
	else
		print_str "WARNING: impossible to add buckling to keep bond_lenght constant"
		print_str "         alat=$alat > alat_0=$ALAT"
	fi
	
	#Cycle over different bucklings
	for buck in $blist; do #0.06 0.08 0.10 0.12 0.14
		let TAB_C++
		set_tab $TAB_C
		POS_2="N       0.500000000   0.288675100   $buck"
		buck_bohr=`echo $buck*$alat | bc -l`
		buck_bohr=`printf %.4f $buck_bohr`

		DESCRIPT="a${alat}_b${buck}"

		print_str "Running scf calculation alat=${alat} (bohr)   buckling=${buck_bohr} (bohr)" "title" $YELLOW

		IN=${PREFIX}_script.${DESCRIPT}_scf.in
		OUT=${PREFIX}_script.${DESCRIPT}_scf.out

		#Set k_point for input
		KPT_MODE="K_POINTS {automatic}"
		KPT_LIST="12 12 1 1 1 1"

		#Print input for scf
		print_in_pw scf
		#Run command for QE
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

		#Extract total energy from output and print it in a two-coloumn file with the buckling
		ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
		echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

		###################################################################################
		#Run band calculation
		echo -e "\n"
		print_str "Running band structure calculation for high-symmetry path of ibrav=$ibrav" "title" $YELLOW
		IN=${PREFIX}_script.${DESCRIPT}_nscf-band.in
		OUT=${PREFIX}_script.${DESCRIPT}_nscf-band.out

		KPT_MODE=""
		KPT_LIST="`cat High_symm/$ibrav.kpt`"

		print_in_pw nscf
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

		#Make the band plot 
		BAND_OUT="bands_b${buck}_plotted.dat"
		do_command "qepp_plotband.x $OUT $BAND_OUT" "null"  $BRIGHT_GREEN
		do_command "gnuplot -e FILE='${BAND_OUT}' -e NBND=$nbnd -e OUTNAME='${DESCRIPT}_bands.pdf' bands.gnu" ""  $BRIGHT_GREEN


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

			IN=${PREFIX}_script.${DESCRIPT}_nscf-kpt${N_KPT_MINGAP}to${app}.in
			OUT=${PREFIX}_script.${DESCRIPT}_nscf-kpt${N_KPT_MINGAP}to${app}_b${buck}.out

			KPT_MODE="K_POINTS {crystal_b}"
			KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

			print_in_pw nscf
			do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

			BAND_OUT="${DESCRIPT}_kpt${N_KPT_MINGAP}to${app}_band.dat"
			do_command "qepp_plotband.x $OUT $BAND_OUT" "null" $BRIGHT_GREEN
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

		IN=${PREFIX}_script.${DESCRIPT}_nscf-opt.in
		OUT=${PREFIX}_script.${DESCRIPT}_nscf-opt.out

		KPT_MODE="K_POINTS {automatic}"
		KPT_LIST="15 15 1 0 0 0"

		print_in_pw nscf
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN


		#Run pw2gw
		IN=${PREFIX}_script.pw2gw.in
		OUT=${PREFIX}_script.${DESCRIPT}_pw2gw.out

		print_in_pw2gw
		do_command "$RUN_COMMAND $BIN_DIR/pw2gw_new.x" "date io" $BRIGHT_GREEN

		do_command "eps_average.x epsX.dat epsY.dat" "null" $BRIGHT_GREEN
		do_command "apply_kk_im.x averaged.dat real_xy.dat 1" "null" $BRIGHT_GREEN

		EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
		VACUUM=`echo "$cdim3 * $alat" | bc -l`
		ALFA0=`echo "$EPS0 * $VACUUM / (4 * $PI)" | bc -l`
		ALFA0=`printf %.5f $ALFA0`

		echo -e "\n${TAB}ALFA0 = $ALFA0"


		####################################################################################
		# Report and final data manipulation
		print_str "Final report for alat=${alat} (bohr)   buckling = $buck_bohr (bohr)" "title" $YELLOW
		print_str "m_h = $MASS_VB    m_e = $MASS_CB   red_mass = $MU" "sub" $CYAN
		print_str "alfa_0 = $ALFA0" "sub" $CYAN

		Aex=`echo "1 / $MU" | bc -l`
		Rex=`echo "$MU / 2" | bc -l`
		#interpolate from Model_2Dexc.plot
		X=`echo "4*$PI*$ALFA0/$Aex" | bc -l`

		X_N=`echo "($X - 0.101) / 0.001 + 1" | bc -l`
		X_N=`printf %.0f $X_N`
		Eb_Rex=`cat model_2Dexc.plot | head -n $X_N | tail -n 1| cut -d$'\t' -f2`
		rex_Aex=`cat model_2Dexc.plot | head -n $X_N | tail -n 1| cut -d$'\t' -f3`

		Eb=`echo "$Eb_Rex * $Rex * $HA_to_EV" | bc -l`
		rex=`echo "$rex_Aex * $Aex" | bc -l`
		print_str "EXC:  Binding = $Eb   radius = $rex" "sub" $CYAN
	done
	let TAB_C--
done
let TAB_C--


echo -e "\nEnd at: $(date)"
echo -e "Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
echo -e "\n\n"




exit









