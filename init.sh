ANSII_RESET="\x1B[0m"
ANSII_F_BLACK="\x1B[30m"
ANSII_F_RED="\x1B[31m"
ANSII_F_GREEN="\x1B[32m"
ANSII_F_YELLOW="\x1B[33m"
ANSII_F_BLUE="\x1B[34m"
ANSII_F_MAGENTA="\x1B[35m"
ANSII_F_CYAN="\x1B[36m"
ANSII_F_WHITE="\x1B[37m"

ANSII_F_BRIGHT_RED="\x1B[38;2;255;0;0m"
ANSII_F_BRIGHT_GREEN="\x1B[38;2;0;255;0m"
ANSII_F_BRIGHT_BLUE="\x1B[38;2;0;0;255m"
ANSII_F_BRIGHT_CYAN="\x1B[38;2;0;255;255m"

export RESET=$ANSII_RESET
export BLACK=$ANSII_F_BLACK
export RED=$ANSII_F_RED
export GREEN=$ANSII_F_GREEN
export YELLOW=$ANSII_F_YELLOW
export BLUE=$ANSII_F_BLUE
export MAGENTA=$ANSII_F_MAGENTA
export CYAN=$ANSII_F_CYAN
export WHITE=$ANSII_F_WHITE

export BRIGHT_RED=$ANSII_F_BRIGHT_RED
export BRIGHT_GREEN=$ANSII_F_BRIGHT_GREEN
export BRIGHT_BLUE=$ANSII_F_BRIGHT_BLUE
export BRIGHT_CYAN=$ANSII_F_BRIGHT_CYAN

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
		echo -e "${TAB}${RED}***${RESET} ${3}${1} ${RED}***${RESET}"
	else
		echo -e "${TAB}${3}${1}${RESET}"
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
			echo -e "${TAB}  ${3}${COMMAND} < ${IN} > ${OUT}${RESET}"
			${COMMAND} < ${IN} > ${OUT}
		else
			echo -e "${TAB}  ${3}${COMMAND} > ${OUT}${RESET}"
			${COMMAND} > ${OUT} 2>&1
		fi
	else
		if [[ `echo $2 | grep -i "null"` != "" ]]; then
			echo -e "${TAB}  ${3}${COMMAND} > /dev/null 2>&1${RESET}"
			$COMMAND  > /dev/null 2>&1
		else
			echo -e "${TAB}  ${3}${COMMAND}${RESET}"
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
