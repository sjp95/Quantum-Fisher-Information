# Quantum Fisher Information (QFI) at Finite Temperature

## Overview
Quantum Fisher Information (QFI) quantifies the ultimate precision bound for estimating a parameter encoded in a quantum state. For thermal (finite-temperature) states, the density matrix is
$$\rho = \frac{e^{-\beta H}}{Z},\quad \beta = \frac{1}{k_B T},\quad Z = \mathrm{Tr}(e^{-\beta H}).$$
In the eigenbasis $\{|n\rangle\}$ of $H(\theta)$ with eigenvalues $\{E_n\}$ and populations $p_n = \frac{e^{-\beta E_n}}{Z}$, the finite-temperature QFI is
$$F_Q = \sum_{n,m} \frac{2(p_n - p_m)^2}{p_n + p_m}\,\big|\langle n|\partial_\theta \rho|m\rangle\big|^2.$$
If the parameter is encoded by a unitary $\rho(\theta)=e^{-i\theta O}\rho e^{i\theta O}$, this simplifies to
$$F_Q = 2\sum_{n,m} \frac{(p_n - p_m)^2}{p_n + p_m}\,\big|\langle n|O|m\rangle\big|^2.$$

For many-body systems, the **QFI density** is $f_Q = \frac{F_Q}{N}$ (with $N$ the number of spins/particles). It serves as an entanglement witness: for $k$-producible states,
$$f_Q \le k,$$
so observing $f_Q > k$ certifies at least $(k+1)$-partite entanglement. For spin-1/2 systems, $f_Q > 1$ implies multipartite entanglement. QFI captures both quantum and thermal fluctuations and is widely used to study metrological sensitivity and critical behavior in many-body systems.

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
