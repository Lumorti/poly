#!/bin/bash

# Variables
cores=8
iters=100000

# User-specified arrays
#prime_numbers=(2 3 4 5 7 8 9 11 13 16 17 19)
#non_prime_numbers=(4 6 8 9 10 12 14 15 18 20)
prime_numbers=(7 8 9 11 13 17 19)
#non_prime_numbers=(4 6 8 9 10 12 14 15 16 18 20)
data_directory="data"

# Function to run the command with given d and n_limit
run_mub_command() {
    local d=$1
    local n_limit=$2

    # Loop for n values
    for ((n=2; n<=n_limit; n++)); do
        ./mub -f -c ${cores} -d ${d} -n ${n} -v 2 -i ${iters} -a 0.9 -b 1e-15 | tee "data/d${d}n${n}.log"
    done
}

# Function to display the line containing "iteration" from a log file
display_iteration_line() {
    local file=$1
    grep "iteration" "$file"
}

# Loop for prime numbers (and their powers)
#for d in "${prime_numbers[@]}"; do
    #n_limit=$((d + 2))
    #run_mub_command $d $n_limit
#done

# Loop for non-prime numbers
#for d in "${non_prime_numbers[@]}"; do
    #n_limit=4
    #run_mub_command $d $n_limit
#done

# Main loop to read all log files in the data directory
for file in "${data_directory}"/*; do
    if [ -f "$file" ]; then
        echo "File: $file"
        display_iteration_line "$file"
        echo
    fi
done

