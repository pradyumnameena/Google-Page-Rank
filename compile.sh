#!/bin/bash
if [ "$1" == "s" ]; then
	clang++ seq.cpp -o mr-pr-seq.o
elif [ "$1" == "c" ]; then
	mpicc mpi_p2p_b.c
elif [ "$1" == "m" ]; then
	mpic++ mr-pr-mpi.cpp -o mr-pr-mpi.o
elif [ "$1" == "mb" ]; then
	clang++ mr-pr-mpi-base.cpp -o mr-pr-mpi-base.o
else
	echo "Supply appropriate and correct arguments"
	exit
fi 
