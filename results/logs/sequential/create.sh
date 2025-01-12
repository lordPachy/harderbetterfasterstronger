#!/bin/bash

printf ".....................................................\n"
printf "SEQUENTIAL LOG FOLDERS CREATION\n"
printf ".....................................................\n"
for seq_length in 10 25 30 40; do
mkdir seq_length_${seq_length}
for num_patterns in 10 15 20 25; do
mkdir seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
mkdir seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
fi
done
done
done
printf "END OF SEQUENTIAL LOG FOLDERS CREATION\n"
printf ".....................................................\n"
