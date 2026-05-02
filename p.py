import glob
import re

N1 = 12
N2 = 100

pattern = f"Data/QFI_{N1}_{N2}_*.dat"
files = glob.glob(pattern)

def extract_last_number(filename):
    # Extract the last number before .dat
    m = re.search(rf"Data/QFI_{N1}_{N2}_(\d+)\.dat", filename)
    if m:
        return int(m.group(1))
    else:
        return None

# Filter and sort files by last number
sorted_files = sorted(
    [(f, extract_last_number(f)) for f in files if extract_last_number(f) is not None],
    key=lambda x: x[1]
)

output_file = f"Data/QFI_{N1}_{N2}_ALL.dat"

with open(output_file, "w") as fout:
    for fname, idx in sorted_files:
        print(f"Adding file: {fname}")

        with open(fname, "r") as fin:
            for line in fin:
                # Skip empty/comment lines if needed
                if line.strip() == "" or line.strip().startswith("#"):
                    continue
                fout.write(line)

print(f"\nMerged file written to: {output_file}")
