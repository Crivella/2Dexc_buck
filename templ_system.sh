#Input file

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
