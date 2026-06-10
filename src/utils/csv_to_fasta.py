import csv

def csv_to_fasta(csv_path, fasta_path, seq_col="designed_sequence", name_col="ID"):
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
csv_to_fasta("final_designs_metrics_10000.csv", "boltzgen_peptides.fasta")