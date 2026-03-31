#!/bin/sh

# Array to store executables and their arguments
commands=(
    "./set_w0_a jack ../../data/jackknife/ ../../data/fit_all_beta ../data_ens_B64"
    "./set_w0_a jack ../../data/jackknife/ ../../data/fit_all_beta ../data_ens_C80"
    "./set_w0_a jack ../../data/jackknife/ ../../data/fit_all_beta ../data_ens_D96"
    "./set_w0_a jack ../../data/jackknife/ ../../data/fit_all_beta ../data_ens_E112"
    "./fit_cont_w0_fpi jack ../../data/jackknife/ ../../data/fit_all_beta"
)

# Array to store failed commands
failed=()

# Function to run a command and capture failures
run_command() {
    local cmd="$1"
    echo "Running: $cmd"
    eval "$cmd"
    if [ $? -ne 0 ]; then
        failed+=("$cmd")
    fi
}

# Run commands sequentially
for cmd in "${commands[@]}"; do
    run_command "$cmd"
done

# Print results
echo "===================="
echo "Results:"
if [ ${#failed[@]} -eq 0 ]; then
    echo "All commands ran successfully."
else
    echo "Failed commands:"
    for cmd in "${failed[@]}"; do
        echo "$cmd"
    done
fi