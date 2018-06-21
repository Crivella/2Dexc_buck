#!/bin/bash

LIST=`grep "\<prefix\>" *.xml`

while read -r line; do 
	file=`echo $line | cut -d: -f1`
	l2=`echo $line | cut -d">" -f2 | cut -d"<" -f1`
	echo $l2
	a=`echo $l2 | cut -d_ -f 2 | tr -dc '[0-9].-'`
	a=`printf %.5f $a`
	b=`echo $l2 | cut -d_ -f 3 | tr -dc '[0-9].-'`
	b=`printf %.5f $b`; pref=`echo $l2 | cut -d_ -f 1`
	sed -i s/$l2/${pref}_a${a}_b${b}/g $file
done <<< "$LIST"


exit










