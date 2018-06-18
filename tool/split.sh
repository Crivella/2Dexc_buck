#!/bin/bash

FILE=`cat $1`

count=0
while read -r line; do
        if [[ `echo $line | grep "#"` != "" ]]; then
                HEADER=`echo -e "$HEADER\n$line"`
        else
                NAME_app=`echo $line | tr -s " " | cut -d " " -f 1`
                NAME="${NAME_app}_$1"
                if [[ `echo $count` == 0 ]]; then
                        echo "" > $NAME
                        echo "$HEADER" >> $NAME
                fi
                echo "      $line" >> $NAME
                let count++
                if [[ `echo $line` == "" ]]; then
                        count=0
                fi
        fi


done <<< "$FILE"


exit








