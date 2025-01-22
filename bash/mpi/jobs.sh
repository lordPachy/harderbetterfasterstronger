#!/bin/bash

printf ".....................................................\n"
printf "MPI JOB CREATION\n"
printf ".....................................................\n"
for num_nodes in 2 4 8; do
mkdir nodes_${num_nodes}
for seq_length in 10 15 20; do
mkdir nodes_${num_nodes}/seq_length_${seq_length}
for num_patterns in 10 15 20 25; do
mkdir nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}
for pattern_mean_length in 4 6 8 10 12 20; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
mkdir nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}
for test_n in {1..10}; do
mkdir nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}
printf "universe = parallel\nexecutable = scripts/mpi/openmpiscript.sh\narguments = scripts/mpi/align_mpi\nshould_transfer_files = YES\ntransfer_input_files = scripts/mpi/align_mpi\nwhen_to_transfer_output = on_exit_or_evict\nlog = logs/mpi/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}\noutput = logs/mpi/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.out.\$(NODE)\nerror = logs/mpi/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.err.\$(NODE)\nmachine_count = ${num_nodes}\nrequest_cpus = 1\ngetenv = True\nqueue" > nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.job
done
fi
done
done
done
done
printf "END OF MPI JOB CREATION\n"
printf ".....................................................\n"
