#!/bin/bash

printf ".....................................................\n"
printf "MPI SCRIPTS CREATION\n"
printf ".....................................................\n"
mkdir ~/scripts/mpi
cp ~/Makefile ~/scripts/mpi
cp ~/rng.c ~/scripts/mpi
cp ~/align_mpi.c ~/scripts/mpi
cp ~/openmpiscript.sh ~/scripts/mpi
cd ~/scripts/mpi
make align_mpi
printf "END OF MPI SCRIPTS CREATION\n"
printf ".....................................................\n"
