#!/bin/bash

FILE=`cat $1`

count=0
while read -r line; do
        if [[ `echo $line | grep "#"` != "" ]]; then
		if [[ `echo $HEADER` != "" ]]; then
                	HEADER=`echo -e "$HEADER\n$line"`
		else
			HEADER=`echo -e "$line"`
		fi
        else
		if [[ $line != "" ]]; then
		        NAME_app=`echo $line | tr -s " " | cut -d " " -f 1`
		        NAME="CD_${NAME_app}_$1"
		        if [[ `echo $count` == 0 ]]; then
		                echo "$HEADER" > $NAME
		        fi
		        echo "      $line" >> $NAME
		        let count++
                else
                        count=0
                fi
        fi


done <<< "$FILE"


exit








