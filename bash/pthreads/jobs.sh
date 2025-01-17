#!/bin/bash

printf ".....................................................\n"
printf "PTHREADS JOB CREATION\n"
printf ".....................................................\n"
for num_threads in 2 4 8 16; do
mkdir threads_${num_threads}
for seq_length in 10 15 20; do
mkdir threads_${num_threads}/seq_length_${seq_length}
for num_patterns in 10 15 20 25; do
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
for test_n in {1..10}; do
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}
printf "universe = vanilla\nlog = logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.log\noutput = logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.out\nerror = logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.err\nrequest_gpus = 0\ngetenv = True\nqueue" > threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.sub
done
fi
done
done
done
done
printf "END OF PTHREADS JOB CREATION\n"
printf ".....................................................\n"
