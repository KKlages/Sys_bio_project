import os

# Mapping from SRR accession to experimental code
srr_to_code = {
    "SRR6480271": "c0n1",
    "SRR6480272": "c1n1",
    "SRR6480273": "c0n0",
    "SRR6480274": "c0n1",
    "SRR6480275": "c1n1",
    "SRR6480276": "c1n0",
    "SRR6480277": "c1n1",
    "SRR6480278": "c0n1",
    "SRR6480279": "c0n0",
    "SRR6480280": "c0n0",
    "SRR6480281": "c1n0",
    "SRR6480282": "c1n0"
}

# Set your results directory here
results_dir = ""

for srr, code in srr_to_code.items():
    quant_dir = os.path.join(results_dir, srr, "kallisto_quant")
    old_file = os.path.join(quant_dir, "abundance.tsv")
    new_file = os.path.join(quant_dir, f"{code}.tsv")
    if os.path.exists(old_file):
        os.rename(old_file, new_file)
        print(f"Renamed {old_file} -> {new_file}")
    else:
        print(f"File not found: {old_file}")

