#!/bin/bash

# Ensure output directory exists
mkdir -p ../Data/qfi

# Set step size (0.05) and range (0.0 to 3.0)
start=0
end=150
step=2

echo "Starting QFI simulations for 18-site Trimer Model..."

for (( hz_int=start; hz_int<=end; hz_int+=step )); do
    # Convert integer back to floating point for python argument
    hz=$(awk -v val="$hz_int" 'BEGIN {printf "%.2f", val / 100.0}')
    
    echo "Running simulation for hz = $hz (File index: qfi_${hz_int}.dat)"
    python3 qfi.py --hz "$hz" --J0 1.0
done

echo "All simulations completed successfully!"