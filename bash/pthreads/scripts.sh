#!/bin/bash

printf ".....................................................\n"
printf "PTHREADS SCRIPTS CREATION\n"
printf ".....................................................\n"
for num_threads in 2 4 8 16; do
mkdir ~/scripts/pthreads/num_threads_${num_threads}
cp ~/Makefile ~/scripts/pthreads/num_threads_${num_threads}
cp ~/rng.c ~/scripts/pthreads/num_threads_${num_threads}
cp ~/align_pthreads.c ~/scripts/pthreads/num_threads_${num_threads}
cd ~/scripts/pthreads/num_threads_${num_threads}
make align_pthread num_threads=${num_threads}
done
printf "END OF PTHREADS SCRIPTS CREATION\n"
printf ".....................................................\n"
