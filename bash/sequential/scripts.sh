#!/bin/bash

printf ".....................................................\n"
printf "SEQUENTIAL SCRIPT CREATION\n"
printf ".....................................................\n"
cp ~/Makefile ~/scripts/sequential
cp ~/rng.c ~/scripts/sequential
cp ~/align.c ~/scripts/sequential
cd ~/scripts/sequential
make align_seq
printf "END OF SEQUENTIAL SCRIPT CREATION\n"
printf ".....................................................\n"
