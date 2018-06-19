#!/bin/bash

export LC_NUMERIC=en_US.UTF-8
PREFIX=$1

if [[ $PREFIX == "" ]]; then
	echo "NO PREFIX SPECIFIED"
	exit
fi

LIST=`find -name "$1*"`

while read -r line; do
PT=`echo $line | rev | cut -d "/" -f 2-30 | rev`
a=`echo $line | cut -d "_" -f 2 | tr -dc '[0-9].-'`
b=`echo $line | cut -d "_" -f 3`
def=`echo $b | cut -d "." -f 3`
if [[ $def == "" ]]; then
	def=`echo $line | cut -d "_" -f 4-20`
	if [[ $def != "" ]]; then
		deflim=_
	else
		deflim=""
	fi
else
	b=`echo $b | cut -d "." -f 1,2`
	deflim=.
fi
b=`echo $b | tr -dc '[0-9].-'`
a=`printf %.5f $a`
b=`printf %.5f $b`

#def=`echo $line | cut -d "_" -f 4-20`

name="${PT}/${PREFIX}_a${a}_b${b}${deflim}$def"

echo "$line --> $name"
mv $line $name
done <<< "$LIST"



exit
