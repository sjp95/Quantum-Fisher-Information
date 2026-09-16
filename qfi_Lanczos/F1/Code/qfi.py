import os

# Set OpenMP environment variables before importing C++ extensions
os.environ["OMP_PROC_BIND"] = "spread"
os.environ["OMP_PLACES"] = "threads"

import argparse
import numpy as np
import qkrylov

def build_trimer_hamiltonian(N, J0, theta, hz):
    """
    Constructs the 18-site Trimer Heisenberg Hamiltonian with PBC and external field hz along Sz.
    J1 = J0 * cos(theta), J2 = J0 * sin(theta)
    Pattern: J1, J1, J2, J1, J1, J2...
    """
    basis = qkrylov.SpinHalfBasis(N, dtype=np.float64)
    site = qkrylov.SpinHalfSite(dtype=np.float64)
    
    J1 = J0 * np.cos(theta)
    J2 = J0 * np.sin(theta)
    
    os_ham = qkrylov.OpSum(dtype=np.float64)
    
    # Exchange interactions with Periodic Boundary Conditions (PBC)
    for i in range(N):
        j = (i + 1) % N  # Connects site N-1 back to site 0
        J_i = J1 if (i % 3 != 2) else J2
        
        os_ham += J_i * 1.0, "Sz", i, "Sz", j
        os_ham += J_i * 0.5, "Sp", i, "Sm", j
        os_ham += J_i * 0.5, "Sm", i, "Sp", j

    # External magnetic field along Sz: -hz * sum_i Sz_i
    for i in range(N):
        os_ham += -hz, "Sz", i

    H = qkrylov.MatrixFreeHamiltonian(basis, site, os_ham, dtype=np.float64)
    return H, basis, site

def compute_qfi_and_magnetization(H, basis, site, N, q=np.pi):
    """
    Computes ground state, QFI density (F_Q / N), and magnetization density (m_z).
    """
    result = qkrylov.lanczos_ground_state(H)
    
    # Extract ground state eigenvector
    if hasattr(result, "eigenvector"):
        psi0 = result.eigenvector
    elif hasattr(result, "vector"):
        psi0 = result.vector
    else:
        psi0 = getattr(result, "psi", result[1] if isinstance(result, tuple) else None)

    # 1. Compute QFI (Sx-Sx correlations at momentum q)
    os_Oq = qkrylov.OpSum(dtype=np.float64)
    for j in range(N):
        phase = np.cos(q * j)  # (-1)^j for q = pi
        os_Oq += phase, "Sx", j

    Oq = qkrylov.MatrixFreeHamiltonian(basis, site, os_Oq, dtype=np.float64)
    psi1 = Oq.apply(psi0)

    exp_Oq = np.vdot(psi0, psi1)
    exp_Oq2 = np.vdot(psi1, psi1)

    qfi_density = 4.0 * np.real(exp_Oq2 - (exp_Oq ** 2)) / N

    # 2. Compute Magnetization density along Sz: m_z = <GS| sum_j Sz_j |GS> / N
    os_Mz = qkrylov.OpSum(dtype=np.float64)
    for j in range(N):
        os_Mz += 1.0, "Sz", j

    Mz_op = qkrylov.MatrixFreeHamiltonian(basis, site, os_Mz, dtype=np.float64)
    psi_mz = Mz_op.apply(psi0)
    mag_density = np.real(np.vdot(psi0, psi_mz)) / N

    return qfi_density, mag_density

def main():
    parser = argparse.ArgumentParser(description="Calculate QFI & Magnetization for 18-site Trimer Model.")
    parser.add_argument("--hz", type=float, required=True, help="External magnetic field along Sz")
    parser.add_argument("--J0", type=float, default=1.0, help="Exchange coupling constant J0")
    parser.add_argument("--theta", type=float, default=np.pi / 4.0, help="Angle theta for trimer coupling")
    args = parser.parse_args()

    N = 18
    q = np.pi

    H, basis, site = build_trimer_hamiltonian(N, args.J0, args.theta, args.hz)
    qfi_val, mag_val = compute_qfi_and_magnetization(H, basis, site, N, q=q)

    # Output directories
    qfi_dir = "../Data/qfi"
    mag_dir = "../Data/Magnetization"
    os.makedirs(qfi_dir, exist_ok=True)
    os.makedirs(mag_dir, exist_ok=True)

    hz_int = int(round(args.hz * 100))
    
    qfi_file = os.path.join(qfi_dir, f"qfi_{hz_int}.dat")
    mag_file = os.path.join(mag_dir, f"mag_{hz_int}.dat")

    # Save QFI data
    with open(qfi_file, "w") as f:
        f.write("# hz\tQFI_density\n")
        f.write(f"{args.hz:.6f}\t{qfi_val:.8f}\n")

    # Save Magnetization data
    with open(mag_file, "w") as f:
        f.write("# hz\tMagnetization_density\n")
        f.write(f"{args.hz:.6f}\t{mag_val:.8f}\n")

    print(f"hz = {args.hz:.4f} (int: {hz_int}) | QFI/N = {qfi_val:.8f} | m_z = {mag_val:.8f}")

if __name__ == "__main__":
    main()