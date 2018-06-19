#!/bin/bash

FILE=`cat $1`
#SYS=`cat system.sh`
source ./system.sh

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
DIST_0=`printf %.5f $DIST_0`

NAME="SD_$1"

count=0
while read -r line; do
        if [[ `echo $line | grep "#"` != "" ]]; then
		if [[ $HEADER != "" ]]; then
                	HEADER=`echo -e "$HEADER\n$line"`
		else
			HEADER=`echo -e "$line"`
		fi
        else
		if [[ $line != "" ]]; then
		        DIST=`echo $line | tr -s " " | cut -d " " -f 3`
			if [[ $DIST == $DIST_0 ]]; then
				if [[ $count == 0 ]]; then
				        echo "$HEADER" > $NAME
				fi
				echo "      $line" >> $NAME
				let count++
			fi
                fi
        fi


done <<< "$FILE"


exit








