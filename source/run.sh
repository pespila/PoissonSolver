#!/bin/bash

# echo
# echo Bash-Script started on `date`
# echo
# make
# ./clean
# dir=/home/michael/Documents/PoissonSolver/tests
# if [ -e $dir ]
# 	then
# 		echo
# 		echo calculation starts
# 	else 
# 		mkdir $dir
# 		echo
# 		echo calculation starts
# fi

# #cp *.txt $dir

# file=dat.txt
# if [ -e $file ]
# 	then
# 		rm $file
# 		echo
# 	else echo
# fi

# for ((i=4;i<=6;i++))
# do
# 	if [ -e dat_$i.txt ]; then
# 		rm dat_$i.txt
# 	fi
# 	for ((j=4096;j<=4096;j=j*2))
# 	do
# 		./poissonSolver $j $i 1 >> dat_$i.txt
# 	done
# done

# mv *.txt $dir
# rm poissonSolver

gnuplot -p plot.gpl