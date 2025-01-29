#!/bin/bash

printf ".....................................................\n"
printf "MPI + PTHREADS LOG FOLDERS CREATION\n"
printf ".....................................................\n"
for num_nodes in 4 8; do
mkdir nodes_${num_nodes}
for num_threads in 2 4 8; do
mkdir nodes_${num_nodes}/threads_${num_threads}
for seq_length in 10 15 20; do
mkdir nodes_${num_nodes}/threads_${num_threads}/seq_length_${seq_length}
for num_patterns in 10 15 20 25; do
mkdir nodes_${num_nodes}/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
mkdir nodes_${num_nodes}/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
for test_n in {1..10}; do
mkdir nodes_${num_nodes}/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}
done
fi
done
done
done
done
done
printf "END OF MPI + PTHREADS LOG FOLDERS CREATION\n"
printf ".....................................................\n"
