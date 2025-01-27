#!/bin/bash

printf ".....................................................\n"
printf "MPI EXP JOB RUNNING\n"
printf ".....................................................\n"
for num_nodes in 2 4 8; do
for seq_length in 10 15 20; do
for num_patterns in 10 15 20 25; do
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
for test_n in {1..5}; do
condor_submit jobs/mpi_exp/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.job
done

printf "Waiting up to 1 minute...\n"
time=0
check="$(condor_q -nobatch | wc -l)"
while [ $time -lt  120 ] && [ $check -gt 9 ]
do
sleep 5s
((time+=5))
check="$(condor_q -nobatch | wc -l)"
done

condor_rm calzona_2046920
fi
done
done
done
done
printf "END OF MPI EXP JOB RUNNING\n"
printf ".....................................................\n"
