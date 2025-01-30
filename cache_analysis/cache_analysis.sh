#!/bin/bash

printf ".....................................................\n"
printf "CACHE ANALYSIS\n"
printf ".....................................................\n"

# Choosing the different settings
settings=("300 0.1 0.3 0.35 100 5 5 300 150 50 150 80 M 609823"
          "10000 0.35 0.2 0.25 0 0 0 10000 9000 9000 50 100 M 4353435"
          "10000 0.35 0.2 0.25 10000 900 900 10000 900 900 500 250 M 4353435"
          "100000 0.35 0.2 0.25 10000 9000 16 10000 9000 16 50000 25000 M 4353435"
          "100000 0.35 0.2 0.25 0 0 0 100000 100000 0 0 0 M 4353435"
)
setting=0

# Compiling all files
cd ..
make align_seq
make align_pthread num_threads=4
make align_mpi
make align_mpi_bcast
make align_mpi_exp
make align_mpipthread num_threads=4
cd cache_analysis

for input in "${settings[@]}"; do
# Creating the setting folder if it does not exist
mkdir sequential/${setting}
mkdir pthreads/${setting}
mkdir mpi/${setting}
mkdir mpi_exp/${setting}
mkdir mpi_bcast/${setting}
mkdir mpi_pthreads/${setting}

# Performing 10 iterations of each test
for test_n in {1..10}; do
sudo perf stat -o sequential/${setting}/${test_n}.txt -B -e cache-references,cache-misses ../align_seq ${input}
sudo perf stat -o pthreads/${setting}/${test_n}.txt -B -e cache-references,cache-misses ../align_pthread ${input}
sudo perf stat -o mpi/${setting}/${test_n}.txt -B -e cache-references,cache-misses mpirun -n 4 ../align_mpi ${input}
sudo perf stat -o mpi_exp/${setting}/${test_n}.txt -B -e cache-references,cache-misses mpirun -n 4 ../align_mpi_exp ${input}
sudo perf stat -o mpi_bcast/${setting}/${test_n}.txt -B -e cache-references,cache-misses mpirun -n 4 ../align_mpi_bcast ${input}
sudo perf stat -o mpi_pthreads/${setting}/${test_n}.txt -B -e cache-references,cache-misses mpirun -n 2 ../align_mpipthread ${input}
done
$((setting++))
done
printf "END OF CACHE ANALYSIS\n"
printf ".....................................................\n"
