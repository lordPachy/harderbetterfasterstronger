#!/bin/bash

printf ".....................................................\n"
printf "CACHE ANALYSIS\n"
printf ".....................................................\n"
settings=("300 0.1 0.3 0.35 100 5 5 300 150 50 150 80 M 609823"
          "1000 0.35 0.2 0.25 0 0 0 20000 10 0 500 0 M 4353435"
          "10000 0.35 0.2 0.25 0 0 0 10000 9000 9000 50 100 M 4353435"
)
counter=0
for 
for input in "${settings[@]}"; do
echo $input
sudo perf stat -o out_$counter.txt -e cache-misses ../align_seq ${input}
$((counter++))
done
printf "END OF CACHE ANALYSIS\n"
printf ".....................................................\n"
