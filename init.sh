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
ANSII_F_ORANGE="\x1B[38;2;255;85;00m"

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
export ORANGE=$ANSII_F_ORANGE

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
			EX=`echo ${COMMAND} | tr " " '\n' | grep ".x" | rev | cut -d/ -f1 | rev | cut -d. -f1`
			PARALLEL=""
			if [[ `echo ${IN} | grep "_scf."` != "" ]]; then
				EX2="scf"
				PARALLEL=${PARALLEL_scf}
			else if [[ `echo ${IN} | grep "_nscf"` != "" ]]; then
				EX2="nscf"
				PARALLEL=${PARALLEL_nscf}
			fi fi
			IN_app=${IN}
			IN="test.in"
			eval "print_in_${EX} ${EX2}"
			IN=${IN_app}

			C1=""
			if [[ -f ${IN} ]]; then
				C1=`diff -q test.in $IN`
			fi
			if [[ `grep "JOB DONE" ${OUT} 2>/dev/null` == "" ]] || [[ ${C1} != "" ]]; then
				if [[ $C1 != "" ]]; then
					mkdir -p old
					if [[ -f ${IN} ]]; then
						mv ${IN} "old/${IN}__${KEEP_DATE}"
					fi
					if [[ -f ${OUT} ]]; then
						mv ${OUT} "old/${OUT}__${KEEP_DATE}"
					fi
				fi
				mv test.in ${IN}
				echo -e "${TAB}  ${3}${COMMAND} ${PARALLEL} < ${IN} > ${OUT}${RESET}"
				${COMMAND} ${PARALLEL} < ${IN} > ${OUT}
			else
				rm test.in
				print_str "Job already performed. Reading from previous output" "" $ORANGE
			fi
		else
			echo -e "${TAB}  ${3}${COMMAND} > ${OUT}${RESET}"
			${COMMAND} > ${OUT} 2>&1
		fi
	else if [[ `echo $2 | grep -i "null"` != "" ]]; then
		echo -e "${TAB}  ${3}${COMMAND} > /dev/null 2>&1${RESET}"
		$COMMAND  > /dev/null 2>&1
	else
		echo -e "${TAB}  ${3}${COMMAND}${RESET}"
		$COMMAND
	fi fi

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
    nat               = $nat
    ntyp              = $ntyp
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
$ATOM_1
$ATOM_2
$ATOM_3
$ATOM_4
$ATOM_5
$ATOM_6

ATOMIC_POSITIONS {$AP}
$pos_1
$pos_2
$pos_3
$pos_4
$pos_5
$pos_6
$pos_7
$pos_8
$pos_9
$pos_10
$pos_11
$pos_12


$KPT_MODE
$KPT_LIST
EOF
}

function print_in_pw2gw()
{
cat > $IN << EOF
&inputpp
  prefix='$PREFIX',
  outdir='$TMP_DIR',
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

function frun_check()
{
	flist=`echo $FRUN_LIST | tr " " "\n"`
	FRUN_CHECK=0
	while read -r line; do
		list1=`echo "$line" | cut -d: -f1 | tr "," "\n"`
		list2=`echo "$line" | cut -d: -f2 | tr "," "\n"`

		C1=""
		C2=""
		while read -r line1; do
			if (( `echo "$line1 == $alat" | bc -l`  )); then
				C1="y"
			fi
		done <<< "$list1"
		while read -r line2; do
			if (( `echo "$line2 == $buck_a" | bc -l`  )); then
				C2="y"
			fi
		done <<< "$list2"

		if [[ $C1 != "" ]] && [[ $C2 != "" ]]; then
			FRUN_CHECK=1
		fi
	done <<< "$flist"
}

function plist_init()
{
	plist=`echo $POINT_LIST | tr " " "\n"`
	while read -r line; do
		list1=`echo "$line" | cut -d: -f1 | tr "," "\n"`

		while read -r line1; do
			if [[ `echo "${ALAT_LIST}" | grep ${line1}` == "" ]]; then
				alist+=" ${line1}"
			fi
		done <<< "$list1"
	done <<< "$plist"
}

export -f print_in_pw
export -f print_in_pw2gw
export -f do_command
export -f print_str
export -f frun_check
export -f plist_init

















