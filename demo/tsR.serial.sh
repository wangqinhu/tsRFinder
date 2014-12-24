#!/bin/bash

# usage
if [ $# != 1 ] ; then
	echo
	echo "Usage: $0 dirname"
	echo " e.g.: $0 conf"
	echo
	exit
fi

dir=$1

for file in `ls $dir/*.conf`; do
	echo "###############################"
	echo "# $file"
	echo "##############################"
	tsRFinder.pl -i no -c $file
done
