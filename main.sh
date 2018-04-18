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
    disk_io = 'minimal'
    wf_collect        = .true.
/
&system
    ibrav             = $ibrav
    celldm(1)         = $alat
    celldm(3)         =  $cdim3
    nat               =  2
    ntyp              = 2
    ecutwfc           =  $cutoff
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

	IN=${PREFIX}_script.scf_b${buck}.in
	OUT=${PREFIX}_script.scf_b${buck}.out

	#Set k_point for input
	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="12 12 1 1 1 1"

	#Print input for scf
	print_in_pw scf

	#Build commmand for pw.x and execute
	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#Extract total energy from output and print it in a two-coloumn file with the buckling
	ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
	echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE

	###################################################################################
	#Run band calculation
	IN=${PREFIX}_script.nscf_b${buck}.in
	OUT=${PREFIX}_script.nscf_b${buck}.out

	KPT_MODE=""
	KPT_LIST="`cat High_symm/$ibrav.kpt`"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#Make the band plot 
	BAND_OUT="bands_b${buck}_plotted.dat"
	#qepp_plotband.x $OUT $BAND_OUT
	gnuplot -e "FILE='$BAND_OUT'" -e "NBND=$nbnd" plot.gnu

	###################################################################################
	#Calcolate effective mass
	#Determine the position of the smallest direct gap
	smallest_gap.x $OUT > tmp_gap.out 2>&1
	N_KPT_MINGAP=`cat tmp_gap.out | grep "Min gap energy:" -A 1 | tail -n 1 | cut -d# -f2  | tr -dc '0-9,-.'`
	rm tmp_gap.out

	AFTER=$N_KPT_MINGAP
	let AFTER++
	BEFORE=$N_KPT_MINGAP
	let BEFORE--

	if [[ $BEFORE == "0" ]]; then
		NUM_KPT=`cat High_symm/$ibrav.kpt | grep -v "#" | wc -l`
		let NUM_KPT--
		let NUM_KPT--
		let NUM_KPT--
		BEFORE=$NUM_KPT
	fi

	echo "$BEFORE $N_KPT_MINGAP $AFTER"

	START=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $N_KPT_MINGAP | tail -n 1 | tr -s " " | cut -d" " -f 2-4`
	END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $AFTER | tail -n 1 | tr -s " " | cut -d" " -f 2-4`

	#Performe calculation with dense k-point near minimum
	IN=${PREFIX}_script.after_b${buck}.in
	OUT=${PREFIX}_script.after_b${buck}.out

	KPT_MODE="K_POINTS {crystal_b}"
	KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	qepp_plotband.x $OUT after_band.dat

	END=`cat High_symm/$ibrav.kpt | tail -n +3 | head -n $BEFORE | tail -n 1 | tr -s " " | cut -d" " -f 2-4`

	IN=${PREFIX}_script.before_b${buck}.in
	OUT=${PREFIX}_script.before_b${buck}.out

	KPT_MODE="K_POINTS {crystal_b}"
	KPT_LIST="`echo -e "2\n$START 100\n$END 1"`"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	qepp_plotband.x $OUT before_band.dat



	####################################################################################
	#Run pw2gw for optical properties
	#Run nscf calculation
	IN=${PREFIX}_script.nscf-opt_b${buck}.in
	OUT=${PREFIX}_script.nscf-opt_b${buck}.out

	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="30 30 1 0 0 0"

	print_in_pw nscf

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`


	#Run pw2gw
	IN=${PREFIX}_script.pw2gw.in
	OUT=${PREFIX}_script.pw2gw_b${buck}.out

	print_in_pw2gw

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw2gw_new.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`

	#eps_average.x epsX.dat epsY.dat
	#apply_kk_im.x averaged.dat real_xy.dat 1

	EPS0=`cat real_xy.dat | grep -v "#" | head -n 1 | tr -s " " | cut -d" " -f3`
	VACUUM=`echo "$cdim3 * $alat" | bc -l`
	ALFA0=`echo "$EPS0 * $VACUUM / (4 * 3.141592)" | bc -l`
	ALFA0=`printf %.5f $ALFA0`

	echo "alfa0 is    $ALFA0"
done






echo "End at: $(date)"




exit









