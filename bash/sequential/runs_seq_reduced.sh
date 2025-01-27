#!/bin/bash

printf ".....................................................\n"
printf "SEQUENTIAL JOB RUNNING\n"
printf ".....................................................\n"

for seq_length in 10; do
for num_patterns in 10 15 20 25; do
for pattern_mean_length in 4; do
if [ ${seq_length} -ge ${pattern_mean_length} ]; then
for test_n in {1..10}; do
condor_submit jobs/sequential/seq_length_${seq_length}/patterns_${num_patterns}/mean_path_length_${pattern_mean_length}/${test_n}/job.sub -append "executable = scripts/sequential/align_seq" "arguments = $(( 2**${seq_length} )) 0.25 0.25 0.25 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( 2**${num_patterns} )) $(( 2**${pattern_mean_length} )) 16 $(( (2**${seq_length})/2 )) $(( (2**${seq_length})/2 )) M 609823" -append 'requirements = (Machine != "node126.di.rm1")'
done

printf "Waiting up to 1 minute...\n"
time=0
check="$(condor_q -nobatch | wc -l)"
while [ $time -lt  60 ] && [ $check -gt 9 ]
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
printf "END OF SEQUENTIAL JOB RUNNING\n"
printf ".....................................................\n"
