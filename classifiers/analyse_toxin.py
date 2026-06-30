import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

TOXINPRED_DIR = "../results/toxinpred"
RESULTS_DIR = "../results"

# Find all toxinpred CSV files
csv_files = sorted(glob.glob(os.path.join(TOXINPRED_DIR, "*_toxinpred.csv")))

labels = []
pct_toxin = []

for csv_path in csv_files:
    df = pd.read_csv(csv_path)
    n_total = len(df)
    n = (df["Prediction"] == "Toxin").sum()
    pct = 100 * n / n_total if n_total else 0
    labels.append(os.path.basename(csv_path).replace("_toxinpred.csv", ""))
    pct_toxin.append(pct)

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(labels, pct_toxin, color="coral", edgecolor="black")
ax.set_ylabel("% Toxin")
ax.set_title("Percentage of toxin predictions per dataset")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
os.makedirs(os.path.join(RESULTS_DIR, "figures"), exist_ok=True)
plt.savefig(os.path.join(RESULTS_DIR, "figures", "toxinpred_barplot.png"), dpi=150, bbox_inches="tight")
plt.close()
