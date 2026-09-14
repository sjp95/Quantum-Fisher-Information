import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def collect_and_plot():
    data_dir = "../Data/qfi"
    combined_file = os.path.join(data_dir, "qfi_vs_hz_all.dat")

    # Match all integer-indexed qfi files
    file_pattern = os.path.join(data_dir, "qfi_[0-9]*.dat")
    files = glob.glob(file_pattern)

    if not files:
        print(f"No data files found in {data_dir}.")
        return

    data = []
    for filepath in files:
        try:
            raw = np.loadtxt(filepath)
            if raw.ndim == 1 and len(raw) >= 2:
                data.append((raw[0], raw[1]))
        except Exception as e:
            print(f"Skipping {filepath}: {e}")

    # Sort data points strictly by hz
    data.sort(key=lambda x: x[0])
    hz_vals = np.array([x[0] for x in data])
    qfi_vals = np.array([x[1] for x in data])

    # Save consolidated dataset
    with open(combined_file, "w") as f:
        f.write("# hz\tQFI\n")
        for hz, qfi in zip(hz_vals, qfi_vals):
            f.write(f"{hz:.6f}\t{qfi:.8f}\n")

    print(f"Collected {len(data)} points into {combined_file}")

    # Plot results
    plt.figure(figsize=(8, 5))
    plt.plot(hz_vals, qfi_vals, "o-", color="navy", linewidth=1.8, markersize=5, label=r"$q = \pi$")
    plt.title(r"Zero-Temperature QFI vs $h_z$ (18-site Trimer Model)", fontsize=13)
    plt.xlabel(r"External Field $h_z$", fontsize=12)
    plt.ylabel(r"Quantum Fisher Information $F_Q$", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=11)
    plt.tight_layout()

    plt.savefig(os.path.join(data_dir, "qfi_vs_hz.png"), dpi=300)
    plt.show()

if __name__ == "__main__":
    collect_and_plot()