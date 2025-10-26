from pathlib import Path
from collections import defaultdict

def rename_duplicate_files():
    base_dir = Path("/user/kilian.klages/u17306/internship_project/mrna_subproject/Afi_data/PRJNA430527/results")
    name_counts = defaultdict(int)

    for run_dir in base_dir.iterdir():
        kallisto_dir = run_dir / "kallisto_quant"
        if not kallisto_dir.is_dir():
            continue

        for tsv_file in kallisto_dir.glob("c*n*.tsv"):
            base_name = tsv_file.stem
            suffix = tsv_file.suffix  # usually ".tsv"

            name_counts[base_name] += 1
            count = name_counts[base_name]

            if count > 1:
                new_name = f"{base_name}_{count}{suffix}"
                new_path = tsv_file.with_name(new_name)

                try:
                    tsv_file.rename(new_path)
                    print(f"Renamed {tsv_file} -> {new_path}")
                except Exception as e:
                    print(f"Failed to rename {tsv_file}: {e}")

if __name__ == "__main__":
    rename_duplicate_files()

