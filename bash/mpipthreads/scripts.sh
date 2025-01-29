#!/bin/bash

printf ".....................................................\n"
printf "MPI + PTHREADS SCRIPTS CREATION\n"
printf ".....................................................\n"
mkdir ~/scripts/mpipthreads
cp ~/openmpiscript.sh ~/scripts/mpipthread
cd ~
make align_mpipthread2 num_threads=2
make align_mpipthread4 num_threads=4
make align_mpipthread8 num_threads=8
printf "END OF MPI + PTHREADS SCRIPTS CREATION\n"
printf ".....................................................\n"
