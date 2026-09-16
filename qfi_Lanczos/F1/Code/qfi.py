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
    basis = qkrylov.SpinHalfBasis(N,dtype=np.float64)
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

    H = qkrylov.MatrixFreeHamiltonian(basis, site, os_ham,dtype=np.float64)
    return H, basis, site

def compute_qfi_zero_temp(H, basis, site, N, q=np.pi):
    """
    Computes zero-temperature QFI for the Sx-Sx correlations at momentum q = pi.
    O_q = sum_j e^{i * q * j} S^x_j = sum_j (-1)^j * S^x_j
    """
    result = qkrylov.lanczos_ground_state(H)
    
    # Access the ground state eigenvector using result.eigenvector (or fallback to result.vector)
    if hasattr(result, "eigenvector"):
        psi0 = result.eigenvector
    elif hasattr(result, "vector"):
        psi0 = result.vector
    else:
        # Fallback if result itself is a tuple or has alternate naming
        psi0 = getattr(result, "psi", result[1] if isinstance(result, tuple) else None)

    os_Oq = qkrylov.OpSum(dtype=np.float64)
    for j in range(N):
        phase = np.cos(q * j)  # (-1)^j for q = pi
        os_Oq += phase, "Sx", j

    Oq = qkrylov.MatrixFreeHamiltonian(basis, site, os_Oq,dtype=np.float64)
    psi1 = Oq.apply(psi0)

    exp_Oq = np.vdot(psi0, psi1)
    exp_Oq2 = np.vdot(psi1, psi1)

    qfi = 4.0 * np.real(exp_Oq2 - (exp_Oq ** 2))
    return qfi/N

def main():
    parser = argparse.ArgumentParser(description="Calculate QFI for 18-site Trimer Model.")
    parser.add_argument("--hz", type=float, required=True, help="External magnetic field along Sz")
    parser.add_argument("--J0", type=float, default=1.0, help="Exchange coupling constant J0")
    parser.add_argument("--theta", type=float, default=np.pi / 4.0, help="Angle theta for trimer coupling")
    args = parser.parse_args()

    N = 18
    q = np.pi

    H, basis, site = build_trimer_hamiltonian(N, args.J0, args.theta, args.hz)
    qfi_val = compute_qfi_zero_temp(H, basis, site, N, q=q)

    output_dir = "../Data/qfi"
    os.makedirs(output_dir, exist_ok=True)

    # Convert hz to integer key (hz * 100) for filename
    hz_int = int(round(args.hz * 100))
    file_path = os.path.join(output_dir, f"qfi_{hz_int}.dat")

    with open(file_path, "w") as f:
        f.write("# hz\tQFI\n")
        f.write(f"{args.hz:.6f}\t{qfi_val:.8f}\n")

    print(f"hz = {args.hz:.4f} (int: {hz_int}) | QFI = {qfi_val:.8f} -> Saved to {file_path}")

if __name__ == "__main__":
    main()