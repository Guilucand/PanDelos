#!/bin/bash

if [ $# -le 1 ]; then
	echo "Not enough arguments! [ilist, prefix] required";
	exit
fi

ilist="$1"
prefix="$2"

mkdir -p $prefix
crtdir=$PWD

cd $prefix
find -name '*.gbk' -size  0 -print0 | xargs -0 rm --
cd $crtdir

while read line
do
echo $line
id=`echo $line | sed s/^.*\_NC\_/NC\_/g`
id=`echo $id | sed s/^.*\_NZ\_/NZ\_/g`
echo "efetch for $id"
if [ ! -f "${prefix}/${id}.gbk" ]; then
	efetch -db nuccore -id $id -format gbwithparts > ${prefix}/${id}.gbk;
fi
done < "$ilist"
python3 gbk2ig.py ${prefix}/ ${prefix}.faa
