#!/bin/bash

printf ".....................................................\n"
printf "PTHREADS JOB CREATION\n"
printf ".....................................................\n"
for num_threads in {2..8}; do
mkdir threads_${num_threads}
for seq_length in 5 10 25 40; do
mkdir threads_${num_threads}/seq_length_${seq_length}
for num_patterns in 5 10 15 20 25; do
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
for pattern_dev in 2 4 8 10; do
mkdir threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/pattern_dev_${pattern_dev}
printf "universe = vanilla\nlog = ~/logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/pattern_dev_${pattern_dev}/job.log\noutput = ~/logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/pattern_dev_${pattern_dev}/job.out\nerror = ~/logs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/pattern_dev_${pattern_dev}/job.err\nrequest_gpus = 0\ngetenv = True\nqueue" > threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/pattern_dev_${pattern_dev}/job.sub
done
done
done
done
done
printf "END OF PTHREADS JOB CREATION\n"
printf ".....................................................\n"
