# Quantum Fisher Information (QFI) at Finite Temperature

## Overview
Quantum Fisher Information (QFI) quantifies the ultimate precision bound for estimating a parameter encoded in a quantum state. For thermal (finite-temperature) states, QFI captures both quantum and thermal fluctuations and is widely used to study metrological sensitivity and critical behavior in many-body systems. For spin-1/2 models, QFI can be computed from the thermal density matrix and related correlation functions, providing a temperature-dependent measure of parameter sensitivity.

## Code Summary
This C++ code computes QFI for a finite-temperature spin-1/2 system. It uses:
- **Intel MKL** for optimized linear algebra routines.
- **OpenMP** for shared-memory parallelism.

The implementation targets high-performance numerical evaluation of QFI at finite temperature. Ensure MKL and OpenMP are available on your system.

## Required Modules (MKL)
```bash
module load compilers/parallel_studio_xe_2018_update3/compiler/default
module load compilers/parallel_studio_xe_2018_update3/mkl/default
```

## Build and Run (MKL)
```bash
export MKL_NUM_THREADS=10
export OMP_NUM_THREADS=10
g++ -std=c++17 -O3 st1.cpp -o h \
	-I${MKLROOT}/include \
	-L${MKLROOT}/lib/intel64 \
	-Wl,--start-group \
		-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread \
	-Wl,--end-group \
	-liomp5 -lpthread -lm -ldl -fopenmp
```
