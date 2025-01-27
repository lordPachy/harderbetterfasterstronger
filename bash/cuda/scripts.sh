#!/bin/bash

printf ".....................................................\n"
printf "CUDA SCRIPT CREATION\n"
printf ".....................................................\n"
cp ~/Makefile ~/scripts/cuda
cp ~/rng.c ~/scripts/cuda
cp ~/align_cuda.cu ~/scripts/cuda
cd ~/scripts/cuda
make align_cuda
printf "END OF CUDA SCRIPT CREATION\n"
printf ".....................................................\n"
