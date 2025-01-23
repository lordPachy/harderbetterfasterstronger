#!/bin/bash

printf ".....................................................\n"
printf "MPI SCRIPTS CREATION\n"
printf ".....................................................\n"
mkdir ~/scripts/mpi
cp ~/openmpiscript.sh ~/scripts/mpi
cd ~
make align_mpi
printf "END OF MPI SCRIPTS CREATION\n"
printf ".....................................................\n"
