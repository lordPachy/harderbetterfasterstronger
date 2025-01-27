#!/bin/bash

printf ".....................................................\n"
printf "MPI BROADCAST SCRIPTS CREATION\n"
printf ".....................................................\n"
mkdir ~/scripts/mpi_broadcast
cp ~/openmpiscript.sh ~/scripts/mpi_broadcast
cd ~
make align_mpi_bcast
printf "END OF MPI BROADCAST SCRIPTS CREATION\n"
printf ".....................................................\n"
