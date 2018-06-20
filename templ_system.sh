#Input file

FIX_DIST="y"			#If != 'n' add a buckling value to each alat to keep a fixed bond lenght set by ALAT_0
				#and the original atom coordinates
ALAT_0=5.72			#Value of alat that defines the constant distance between atoms
				#The constant distance is determined at run_time from ALAT_0, POS_1 and POS_2
CDIM3_0=6			#Value of celldm(3) that defines the constant vacuum lenght
				#Constant lenght for vacuum of CDIM3_0 * $ALAT_0
ALAT_LIST="5.70 5.72 5.74"	#List of alat to cycle on
BUCKLING_LIST="0.00"		#List of bucklings to cycle on
				#An additional buckling for every alat will be added to also give the result for constant distance

#QE input cards values
prefix=AlN
cutoff=30 	#100
nat=2
ntyp=2
nbnd=20
ibrav=4

KPT_LIST_scf="12 12 1 1 1 1"

DO_NSCF_OPT="y"			#If != 'n' run the nscf calculation for the optical properties
DO_PW2GW="y"			#If != 'n' run the pw2gw calculation for the optical properties
DO_PW2GW_FRUN="n"		#If != 'n' run the pw2gw calculation for the 'force run' list cases
FRUN_LIST=""			#list of nscf/pw2gw calculation to run anyway if encountered in the normal cycle
				#example "4.6:0.1 4.7:0.2 ..." (celldim:buck) separate them by using spaces
				#example "4.6:0.1,0.2,0.3" use , to make a list of cdim/buck
KPT_LIST_pw2gw="102 102 1 0 0 0"

PARALLEL_scf="-npools 2" #-npools 17
PARALLEL_nscf="-npools 2" #-npools 34

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
