#!/bin/bash
export LC_NUMERIC=en_US.UTF-8

FILE=`cat $1`
source ./system.sh


#rm BK_*
while read -r line; do
	if [[ `echo $line | grep "#"` != "" ]]; then
		HEADER=`echo -e "$HEADER\n$line"`
	else
		if [[ $line != "" ]]; then
			#CD=`echo $line | tr -s " " | cut -d " " -f 1`
			buck_app=`echo $line | tr -s " " | cut -d " " -f 2`
			BUCK=`echo $buck_app/$ALAT_0 | bc -l`
			BUCK=`printf %.4f $BUCK`
			NAME="BK_${BUCK}_$1"
			if test -e $NAME; then
				if [[ `cat $NAME` == "" ]]; then
					echo "$HEADER" > $NAME
				fi
			else
				echo "$HEADER" > $NAME
			fi
			echo "      $line" >> $NAME 
		fi
	fi
done <<< "$FILE"


exit











