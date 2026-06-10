import pandas as pd
import os
from transformers import AutoModelForCausalLM
from evaluate import load
import numpy as np

fasta_file = "data/fastas/dit_amp.fasta"

# Load perplexity metric
perplexity = load("perplexity", module_type="metric")

# Read the CSV file
from Bio import SeqIO

fasta_sequences = list(SeqIO.parse(fasta_file, "fasta"))
df = pd.DataFrame({
    'sequence': [str(record.seq) for record in fasta_sequences]
})

sequences = df['sequence'].tolist()

results = perplexity.compute(predictions=sequences, model_id='hugohrban/progen2-medium')
df['progen2_ppl'] = results['perplexities']

df.to_csv(f"{os.path.basename(fasta_file)}_with_progen2_ppl.csv", index=False)

mean_ppl = np.mean(df['progen2_ppl'])
std_ppl = np.std(df['progen2_ppl'])
print(f"Mean PPL for progen2: {mean_ppl:.2f} ± {std_ppl:.2f}")

