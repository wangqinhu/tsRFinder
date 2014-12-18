#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -N tsRFinder
#$ -V

if [ $# != 1 ] ; then
	echo
	echo "Usage: $0 config_file"
	echo " e.g.: $0 tsR.conf"
	echo
	exit
fi

echo -n "Running on: "
hostname
echo "SGE job id: $JOB_ID"
echo -n "Begin time: "
date

echo
echo "||||||||||||||||||||||||||||||||||||||||"
echo

tsRFinder.pl -o yes -c $1

echo
echo "||||||||||||||||||||||||||||||||||||||||"
echo

echo -n "End time: "
date
