csv_paths = ["KPC3_designs/final_designs_metrics_100.csv", "NDM5_designs/final_designs_metrics_100.csv"]
output_path = "all_designs.fasta"

import csv

def csv_to_fasta(csv_path, fasta_path, seq_col="designed_sequence", name_col="id"):
    with open(csv_path, newline='') as csvfile, open(fasta_path, "w") as fasta:
        reader = csv.DictReader(csvfile)
        for row in reader:
            seq = row.get(seq_col, "").strip()
            name = row.get("final_rank", "").strip()
            if seq:
                if not name:
                    name = "unranked"
                fasta.write(f">seq_rank_{name}\n{seq}\n")

# Example usage:
csv_to_fasta(csv_paths[0], "KPC3_designs.fasta")
csv_to_fasta(csv_paths[1], "NDM5_designs.fasta")
