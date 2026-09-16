#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

mkdir -p ../Data/qfi

start=0
end=200
step=2

echo "Starting QFI simulations for 18-site Trimer Model..."

for (( hz_int=start; hz_int<=end; hz_int+=step )); do
    hz=$(awk -v val="$hz_int" 'BEGIN {printf "%.2f", val / 100.0}')
    
    echo "Running simulation for hz = $hz (File index: qfi_${hz_int}.dat)"
    
    # Pass OMP flags inline right when launching Python
    OMP_PROC_BIND=spread OMP_PLACES=threads python3 qfi.py --hz "$hz" --J0 1.0
done

echo "All simulations completed successfully!"