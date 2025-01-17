#!/bin/bash

printf ".....................................................\n"
printf "SEQUENTIAL JOB CREATION\n"
printf ".....................................................\n"
for seq_length in 10 15 20; do
mkdir seq_length_${seq_length}
for num_patterns in 10 15 20 25; do
mkdir seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
mkdir seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
for test_n in {1..10}; do
mkdir seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}
printf "universe = vanilla\nlog = logs/sequential/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.log\noutput = logs/sequential/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.out\nerror = logs/sequential/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.err\nrequest_gpus = 0\ngetenv = True\nqueue" > seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.sub
done
fi
done
done
done
printf "END OF SEQUENTIAL JOB CREATION\n"
printf ".....................................................\n"
