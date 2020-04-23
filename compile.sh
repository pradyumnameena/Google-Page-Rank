#!/bin/bash
if [ "$1" == "s" ]; then
	clang++ seq.cpp -o mr-pr-seq.o
elif [ "$1" == "c" ]; then
	# clang++ mr-pr-cpp.cpp -o mr-pr-cpp.o
	mpicc mpi_p2p_b.c
elif [ "$1" == "m" ]; then
	mpicc mpi_p2p_n.c
elif [ "$1" == "mb" ]; then
	mpicc mpi_c.c
else
	echo "Supply appropriate and correct arguments"
	exit
fi 
