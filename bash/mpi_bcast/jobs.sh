#!/bin/bash

printf ".....................................................\n"
printf "MPI BROADCAST JOB CREATION\n"
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
printf "universe = parallel\n\nexecutable = scripts/mpi_bcast/openmpiscript.sh\n\narguments = align_mpi_bcast $(( 2**${seq_length} )) 0.25 0.25 0.25 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( (2**${seq_length})/2 )) $(( (2**${seq_length})/4 )) M 609823\n\nshould_transfer_files = YES\n\ntransfer_input_files = align_mpi_bcast\n\nwhen_to_transfer_output = on_exit_or_evict\n\noutput = logs/mpi_bcast/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.out.\$(NODE)\nerror = logs/mpi_bcast/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.err.\$(NODE)\nlog = logs/mpi_bcast/nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.log\n\nmachine_count = ${num_nodes}\nrequest_cpus = 1\n\ngetenv = True\nrequirements = (Machine != \"node126.di.rm1\")\nqueue\n" > nodes_${num_nodes}/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.job
done
fi
done
done
done
done
printf "END OF MPI BROADCAST JOB CREATION\n"
printf ".....................................................\n"
