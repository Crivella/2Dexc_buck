#!/bin/bash

source ./ENVIRONMENT_VARIABLES
echo "BIN_DIR:" $BIN_DIR
echo "PSEUDO_DIR:" $PSEUDO_DIR
echo "TMP_DIR:" $TMP_DIR
echo "Parallel command:" $RUN_COMMAND
echo "Started at: " `date`

function print_in()
{
cat > $IN << EOF
&control
    calculation       = '$1'
    restart_mode      = 'from_scratch'
    prefix            = 'AlN'
    tstress           = .true.
    tprnfor           = .true.
    pseudo_dir        = '$PSEUDO_DIR'
    outdir            = '$TMP_DIR'
    verbosity         = 'high'
    disk_io = 'minimal'
    wf_collect        = .true.
/
&system
    ibrav             = 4
    celldm(1)         = $alat
    celldm(3)         =  6
    nat               =  2
    ntyp              = 2
    ecutwfc           =  100
    nbnd              = 20
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
Al  26.982  Al.pz-vbc.upf
N   14.067  N.pz-vbc.upf

ATOMIC_POSITIONS {alat}
Al      0.000000000   0.577350300   0.00
N       0.500000000   0.288675100   $buck

$KPT_MODE
$KPT_LIST
EOF
}

export print_in()

SAVE=Etot_vs_buck.dat
echo -e "# buckling(bohr) Etot(Ry)" > $SAVE

for buck in 0.06 0.08 0.10 0.12 0.14; do
	IN=AlN_script.scf_b${buck}.in
	OUT=AlN_script.scf_b${buck}.out

	alat=5.72

	buck_bohr=`echo $buck*$alat | bc -l`
	buck_bohr=`printf %.4f $buck_bohr`

	KPT_MODE="K_POINTS {automatic}"
	KPT_LIST="12 12 1 1 1 1"

	print_in scf

	echo -e "\tStart: " `date`
	COMMAND="  $RUN_COMMAND $BIN_DIR/pw.x"
	echo -e "\t\t$COMMAND < $IN > $OUT"
	#$COMMAND < $IN > $OUT
	echo -e "\tEnd: " `date`


	ENERGY=`cat $OUT | grep ! | tr -dc '0-9,-.'`
	echo -e "$buck_bohr\t\t$ENERGY" >> $SAVE


	IN=AlN_script.nscf_b${buck}.in
	OUT=AlN_script.nscf_b${buck}.out
done






echo "End at: $(date)"




exit









