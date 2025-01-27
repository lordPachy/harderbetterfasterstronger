#!/bin/bash

printf ".....................................................\n"
printf "MPI EXP SCRIPTS CREATION\n"
printf ".....................................................\n"
mkdir ~/scripts/mpi_exp
cp ~/openmpiscript.sh ~/scripts/mpi_exp
cd ~
make align_mpi_exp
printf "END OF MPI EXP SCRIPTS CREATION\n"
printf ".....................................................\n"
