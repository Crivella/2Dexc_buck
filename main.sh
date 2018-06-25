#!/bin/bash

source ./ENVIRONMENT_VARIABLES

KEEP_DATE=`date +%Y-%m-%d_%H.%M.%S`
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND


source ./init.sh	#File containing general variable and function definitions
			#do_command() print_str() print_in_pw() print_in_pw2gw()
source ./constants.sh	#File containing physical constants definitions
source ./system.sh	#File containing data of the system to analyze


START_TIME=`date +%s.%N`
echo "Started at: " `date`

TAB_C=0
set_tab $TAB_C

SAVE="${prefix}_SAVE.dat"
save_file ${SAVE}
printf "#%10s%12s%11s%13s%10s%10s%10s%12s%9s%9s%7s%13s%14s\n" "alat(bohr)" "buck(bohr)" "dist(bohr)" "Etot(Ry)" "m_h1(au)" "m_h2(au)" "m_e(au)" "m_red(a.u.)" "E_vb(eV)" "E_cb(eV)" "alfa_0" "Exc_b_en(eV)" "Exc_rad(a.u.)" > $SAVE


#Cycle over different celldimension
A1_a=`echo $POS_1 | tr -s " " | cut -d" " -f 1`
A1_x=`echo $POS_1 | tr -s " " | cut -d" " -f 2`
A1_y=`echo $POS_1 | tr -s " " | cut -d" " -f 3`
A1_z=`echo $POS_1 | tr -s " " | cut -d" " -f 4`
A2_a=`echo $POS_2 | tr -s " " | cut -d" " -f 1`
A2_x=`echo $POS_2 | tr -s " " | cut -d" " -f 2`
A2_y=`echo $POS_2 | tr -s " " | cut -d" " -f 3`
A2_z=`echo $POS_2 | tr -s " " | cut -d" " -f 4`
Dx=`echo "$A1_x - $A2_x" | bc -l`
Dy=`echo "$A1_y - $A2_y" | bc -l`
Dz=`echo "$A1_z - $A2_z" | bc -l`
DIST_0=`echo "${ALAT_0} * sqrt($Dx^2 + $Dy^2 + ($Dz * $CDIM3_0)^2)" | bc -l`

alist=${ALAT_LIST}
#Add specific point to cicle through
plist_init

let TAB_C++
for alat_a in ${alist}; do
	set_tab $TAB_C
	alat=`printf %.5f $alat_a`
	cdim3=`echo "${ALAT_0}*${CDIM3_0} / $alat" | bc -l`
	cdim3=`printf %.4f $cdim3`
	print_str "Running alat=${alat}   cdim3=${cdim3}" "title" $YELLOW

	#Introduce buckling value to calculate that gives atom distance equal to DIST_0
	#Fixed atom distance calculation
	
	blist="$BUCKLING_LIST "
	if [[ "$FIX_DIST" != "n" ]]; then
		if (( `echo "$alat <= ${ALAT_0}" | bc -l` )); then
			CHECK=0
			for b in $BUCKLING_LIST; do
				if (( `echo "$Dz >= 0" | bc -l` )); then
					BCK=1
					DIST=`echo "$alat * sqrt($Dx^2 + $Dy^2 + ($Dz + $b)^2 * ($cdim3)^2 )" | bc -l`
				else
					BCK=2
					DIST=`echo "$alat * sqrt($Dx^2 + $Dy^2 + ($Dz - $b)^2 * ($cdim3)^2 )" | bc -l`
				fi
				if (( `echo "sqrt(($DIST - $DIST_0)^2)<0.00005" | bc -l` )); then
					CHECK=1
					break
				fi
			done
			added=""
			if [[ $CHECK == "0" ]]; then
				Arat=`echo "$ALAT_0 / $alat" | bc -l`
				b=`echo "sqrt(($Dx^2 + $Dy^2)*($Arat^2 - 1) + $Arat^2*$Dz) - $Dz" | bc -l`
				b=`printf %0.5f $b`
				added="$b"
				blist+=" $b "
				print_str "Added buckling value $b to loop to make constant distant calculation"
			fi
		else
			print_str "WARNING: impossible to add buckling to keep bond_lenght constant"
			print_str "         alat=$alat > alat_0=${ALAT_0}"
		fi
	fi

	CP=`echo "${plist}" | grep ${alat_a}`
	if [[ $CP != "" ]]; then
		if [[ `echo ${ALAT_LIST} | grep ${alat_a}` == "" ]]; then
			blist=""
		fi
		blist+=`echo "${CP}" | cut -d ":" -f 2 | tr "," " "`
	fi

	#Sort buckling list
	blist=`echo $blist | tr " " "\n" | sort | tr "\n" " "`
	
	#Cycle over different bucklings
	let TAB_C++
	for buck in $blist; do
		set_tab $TAB_C

		buck_a=`printf %.5f $buck`
		if [[ $buck != $added ]]; then
			buck=`echo "$buck * $ALAT_0/$alat" | bc -l`
		fi

		if [[ $BCK=="1" ]]; then
			Z=`echo "$A1_z + $buck" | bc -l`
			Z=`printf %.8f $Z`
			pos_1="$A1_a $A1_x $A1_y $Z"
			pos_2="$A2_a $A2_x $A2_y $A2_z"
		fi
		if [[ $BCK=="2" ]]; then
			Z=`echo "$A2_z + $buck" | bc -l`
			Z=`printf %.8f $Z`
			pos_1="$A1_a $A1_x $A1_y $A1_z"
			pos_2="$A2_a $A2_x $A2_y $Z"
		fi
		buck_bohr=`echo $buck*$alat | bc -l`
		buck_bohr=`printf %.4f $buck_bohr`
		DIST=`echo "sqrt($Dx^2 + $Dy^2 + (sqrt($Dz^2) + $buck )^2 ) * $alat" | bc -l`

		PREFIX="${prefix}_a${alat}_b${buck_a}"

		#SAVE="${PREFIX}_SAVE.dat"
		#printf "#%14s%14s%14s%14s" "alat(bohr)" "buck(bohr)" "dist(bohr)" "Etot(Ry)" > $SAVE

		print_str "Running scf calculation alat=${alat} (bohr)   buckling=${buck_bohr} (bohr)" "title" $YELLOW

		IN=${PREFIX}_scf.in
		OUT=${PREFIX}_scf.out

		#Set k_point for input
		KPT_MODE="K_POINTS {automatic}"
		KPT_LIST=${KPT_LIST_scf}

		#Print input for scf
		#print_in_pw scf
		#Run command for QE
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

		cp -t tmp $IN $OUT

		#Extract total energy from output and print it in a two-coloumn file with the buckling
		ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
		#echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

		###################################################################################
		#Run band calculation
		echo -e "\n"
		print_str "Running band structure calculation for high-symmetry path of ibrav=$ibrav" "title" $YELLOW
		IN=${PREFIX}_nscf-band.in
		OUT=${PREFIX}_nscf-band.out

		KPT_MODE=""
		KPT_LIST="`cat High_symm/$ibrav.kpt`"

		#print_in_pw nscf
		do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

		#Make the band plot 
		BAND_OUT="${PREFIX}_plotted.dat"
		if [[ ! -f $BAND_OUT ]]; then
			do_command "tool/bin/qepp_plotband.x $OUT $BAND_OUT" "null" $BRIGHT_GREEN
		fi
		do_command "gnuplot -e FILE='${BAND_OUT}' -e NBND=$nbnd -e OUTNAME='${BAND_OUT:0: -4}.pdf' bands.gnu" ""  $BRIGHT_GREEN


		###################################################################################
		#Calcolate effective mass
		#Determine the position of the smallest direct gap and valence and conduction band number
		echo -e "\n"
		print_str "Starting the calcuations for the effective mass" $YELLOW
		IN=""
		APP="$OUT"
		OUT="tmp_gap.out"
		do_command "tool/bin/smallest_gap.x $APP" "io null" $BRIGHT_GREEN

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
		MASS_VB_HEAVY=0
		MASS_VB_LIGHT=0
		MASS_CB_HEAVY=0
		MASS_CB_LIGHT=0
		for app in $AFTER $BEFORE; do
			MASS_VB_LIGHT_APP="0"
			MASS_CB_LIGHT_APP="0"
			END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $app | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
			echo -e "${TAB}Running direction $START -> $END"

			IN=${PREFIX}_nscf-kpt${N_KPT_MINGAP}to${app}.in
			OUT=${PREFIX}_nscf-kpt${N_KPT_MINGAP}to${app}.out

			KPT_MODE="K_POINTS {crystal_b}"
			KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

			#print_in_pw nscf
			do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN

			BAND_OUT="${PREFIX}_kpt${N_KPT_MINGAP}to${app}_band.dat"
			if [[ ! -f $BAND_OUT ]]; then
				do_command "tool/bin/qepp_plotband.x $OUT $BAND_OUT" "null" $BRIGHT_GREEN
			fi
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($VB+1)) -e OUTNAME='${BAND_OUT:0: -4}_vb.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
			E_BAND_0_VB=`cat app1.dat`
			MASS_VB_HEAVY_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
			do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+1)) -e OUTNAME='${BAND_OUT:0: -4}_cb.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
			E_BAND_0_CB=`cat app1.dat`
			MASS_CB_HEAVY_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
			rm app.dat
			rm app1.dat

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
				MASS_VB_LIGHT_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
			else
				print_str "  No degeneracy for the valence band"
				print_str "  vb and vb-1 are split more than 1meV"
			fi

			if [[ $DEGEN == "1" ]]; then
				print_str "m_vb = $MASS_VB_HEAVY_APP    m_vb-1 = $MASS_VB_LIGHT_APP" "" $CYAN
				MASS_VB=`echo "$MASS_VB + sqrt(($MASS_VB_HEAVY_APP + $MASS_VB_LIGHT_APP)^2)/4" | bc -l`
				MASS_VB_HEAVY=`echo "$MASS_VB_HEAVY + sqrt(($MASS_VB_HEAVY_APP)^2)/2" | bc -l`
				MASS_VB_LIGHT=`echo "$MASS_VB_LIGHT + sqrt(($MASS_VB_LIGHT_APP)^2)/2" | bc -l`
			else
				print_str "m_vb = $MASS_VB_HEAVY_APP" "" $CYAN
				MASS_VB=`echo "$MASS_VB + sqrt(($MASS_VB_HEAVY_APP)^2)/2" | bc -l`
				MASS_VB_HEAVY=`echo "$MASS_VB_HEAVY + sqrt(($MASS_VB_HEAVY_APP)^2)/2" | bc -l`
			fi

			E_CB=`echo $MIN_LINE | cut -d" " -f $CB`
			E_CB_1=`echo $MIN_LINE | cut -d" " -f $(($CB+1))`
			DIFF=`echo "$E_CB_1 - $E_CB" | bc -l`

			DEGEN=0
			print_str "Looking for degeneracy of the conduction band" "sub"
			print_str "  e(cb) = $E_CB,    e(cb+1) = $E_CB_1  ->  $DIFF"
			if (( `echo "$DIFF<0.001" | bc -l` )); then
				DEGEN=1
				print_str "  Degeneracy found within the 1meV limit.."
				print_str "  running fit of the cb+1 band"
				do_command "gnuplot -e FILE='$BAND_OUT' -e NBND=$(($CB+2)) -e OUTNAME='${BAND_OUT:0: -4}_cb-1.pdf' emass_fit.gnu" "" $BRIGHT_GREEN
				MASS_CB_LIGHT_APP=`echo "(2* $(cat app.dat) / $HA_to_EV * ($alat /(2*$PI))^2)^(-1)" | bc -l`
			else
				print_str "  No degeneracy for the conduction band"
				print_str "  cb and cb+1 are split more than 1meV"
			fi

			if [[ $DEGEN == "1" ]]; then
				print_str "m_cb = $MASS_CB_HEAVY_APP    m_cb+1 = $MASS_CB_LIGHT_APP" "" $CYAN
				MASS_CB=`echo "$MASS_CB + sqrt(($MASS_CB_HEAVY_APP + $MASS_CB_LIGHT_APP)^2)/4" | bc -l`
				MASS_CB_HEAVY=`echo "$MASS_CB_HEAVY + sqrt(($MASS_CB_HEAVY_APP)^2)/2" | bc -l`
				MASS_CB_LIGHT=`echo "$MASS_CB_LIGHT + sqrt(($MASS_CB_LIGHT_APP)^2)/2" | bc -l`
			else
				print_str "m_cb = $MASS_CB_HEAVY_APP" "" $CYAN
				MASS_CB=`echo "$MASS_CB + sqrt(($MASS_CB_HEAVY_APP)^2)/2" | bc -l`
				MASS_CB_HEAVY=`echo "$MASS_CB_HEAVY + sqrt(($MASS_CB_HEAVY_APP)^2)/2" | bc -l`
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

		IN=${PREFIX}_nscf-opt.in
		OUT=${PREFIX}_nscf-opt.out

		KPT_MODE="K_POINTS {automatic}"
		KPT_LIST=${KPT_LIST_pw2gw}

		frun_check

		#print_in_pw nscf
		if [[ "$DO_NSCF_OPT" != "n" ]] || [[  $FRUN_CHECK == 1 ]]; then
			do_command "$RUN_COMMAND $BIN_DIR/pw.x" "date io" $BRIGHT_GREEN
		fi


		#Run pw2gw
		IN=${PREFIX}_pw2gw.in
		OUT=${PREFIX}_pw2gw.out

		#print_in_pw2gw
		if [[ "$DO_PW2GW" != "n" ]] || [[  $FRUN_CHECK == 1 && $DO_PW2GW_FRUN != "n" ]]; then
			do_command "$RUN_COMMAND $BIN_DIR/pw2gw.x" "date io" $BRIGHT_GREEN
		fi

		if [[ ! -d ${PREFIX} ]]; then
			mkdir -p ${PREFIX}
			mv -t ${PREFIX} matrixelements k.dat eps*
		else if [[ -f "matrixelements" ]]; then
			mv -t ${PREFIX} matrixelements k.dat eps*
		fi fi
		cd ${PREFIX}

		if [[ ! -f "averaged.dat" ]]; then
			if [[ -f "epsX.dat" ]]; then
				do_command "../tool/bin/eps_average.x epsX.dat epsY.dat" "null" $BRIGHT_GREEN
			fi
		fi
		if [[ ! -f "real_xy.dat" ]]; then
			if [[ -f "averaged.dat" ]]; then
				do_command "../tool/bin/apply_kk_im.x averaged.dat real_xy.dat 1" "null" $BRIGHT_GREEN
			fi
		fi		

		EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
		VACUUM=`echo "$cdim3 * $alat" | bc -l`
		ALFA0=`echo "($EPS0 - 1) * $VACUUM / (4 * $PI)" | bc -l`
		ALFA0=`printf %.5f $ALFA0`

		print_str "Alfa_0 = ${ALFA0}" "sub" $CYAN
		#echo -e "\n${TAB}ALFA0 = $ALFA0"

		cd ..


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


		printf "%11.5f%12.5f%11.5f%13.7f%10.5f%10.5f%10.5f%12.5f%9.4f%9.4f%7.4f%13.6f%14.6f\n" "$alat" "$buck_bohr" "$DIST" "$ENERGY" "$MASS_VB_HEAVY" "$MASS_VB_LIGHT" "$MASS_CB" "$MU" "$E_BAND_0_VB" "$E_BAND_0_CB" "$ALFA0" "$Eb" "$rex" >> $SAVE
	done
	let TAB_C--
	printf "\n" >> $SAVE
done
let TAB_C--


echo -e "\nEnd at: $(date)"
echo -e "Total run time: " `echo "$(date +%s.%N) - $START_TIME" | bc -l` "s"
echo -e "\n\n"




exit









