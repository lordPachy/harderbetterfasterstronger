#!/bin/bash

printf ".....................................................\n"
printf "PTHREADS JOB RUNNING\n"
printf ".....................................................\n"
for num_threads in 2 4 8 16; do
for seq_length in 10 15 20; do
for num_patterns in 10 15 20 25; do
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
for test_n in {1..10}; do
condor_submit jobs/pthreads/threads_${num_threads}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.sub -append "executable = scripts/pthreads/num_threads_${num_threads}/align_pthread" "arguments = $(( 2**${seq_length} )) 0.25 0.25 0.25 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( (2**${seq_length})/2 )) $(( (2**${seq_length})/2 )) M 609823" -append 'requirements = (Machine != "node126.di.rm1")'
done
fi
done
done
done
done
printf "END OF PTHREADS JOB RUNNING\n"
printf ".....................................................\n"
