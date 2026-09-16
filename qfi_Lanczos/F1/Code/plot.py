import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def load_and_sort_data(data_dir, prefix):
    file_pattern = os.path.join(data_dir, f"{prefix}_[0-9]*.dat")
    files = glob.glob(file_pattern)

    data = []
    for filepath in files:
        try:
            raw = np.loadtxt(filepath)
            if raw.ndim == 1 and len(raw) >= 2:
                data.append((raw[0], raw[1]))
        except Exception as e:
            print(f"Skipping {filepath}: {e}")

    data.sort(key=lambda x: x[0])
    return data

def main():
    qfi_dir = "../Data/qfi"
    mag_dir = "../Data/Magnetization"

    qfi_data = load_and_sort_data(qfi_dir, "qfi")
    mag_data = load_and_sort_data(mag_dir, "mag")

    if not qfi_data or not mag_data:
        print("Data files missing in ../Data/qfi or ../Data/Magnetization.")
        return

    # Consolidate QFI dataset
    qfi_all_file = os.path.join(qfi_dir, "qfi_vs_hz_all.dat")
    with open(qfi_all_file, "w") as f:
        f.write("# hz\tQFI_density\n")
        for hz, val in qfi_data:
            f.write(f"{hz:.6f}\t{val:.8f}\n")

    # Consolidate Magnetization dataset
    mag_all_file = os.path.join(mag_dir, "mag_vs_hz_all.dat")
    with open(mag_all_file, "w") as f:
        f.write("# hz\tMagnetization_density\n")
        for hz, val in mag_data:
            f.write(f"{hz:.6f}\t{val:.8f}\n")

    print(f"Saved combined files:\n  - {qfi_all_file}\n  - {mag_all_file}")

    # Extract arrays for plotting
    hz_qfi, qfi_vals = zip(*qfi_data)
    hz_mag, mag_vals = zip(*mag_data)

    # Plot two subplots (QFI and Magnetization vs hz)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # Subplot 1: QFI
    ax1.plot(hz_qfi, qfi_vals, "o-", color="navy", linewidth=1.5, markersize=4, label=r"$F_Q/N$ ($q=\pi$)")
    ax1.set_ylabel(r"QFI Density ($F_Q/N$)", fontsize=12)
    ax1.set_title(r"18-site Trimer Model ($J_1-J_1-J_2$)", fontsize=14)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(fontsize=11)

    # Subplot 2: Magnetization
    ax2.plot(hz_mag, mag_vals, "s-", color="firebrick", linewidth=1.5, markersize=4, label=r"$m_z$")
    ax2.set_xlabel(r"External Magnetic Field ($h_z$)", fontsize=12)
    ax2.set_ylabel(r"Magnetization Density ($m_z$)", fontsize=12)
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(fontsize=11)

    plt.tight_layout()
    qfi_plot = "../Plot/qfi"
    plt.savefig(os.path.join(qfi_plot, "qfi_mag_vs_hz.pdf"), dpi=300)
    plt.show()

if __name__ == "__main__":
    main()