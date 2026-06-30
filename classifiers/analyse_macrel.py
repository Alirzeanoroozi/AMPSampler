import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_DIR = "../results"

# Find all *_macrel dirs
macrel_dirs = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_macrel")))

labels = []
pct_amp = []

for d in macrel_dirs:
    csv_path = os.path.join(d, "macrel.csv")
    if not os.path.exists(csv_path):
        continue
    df = pd.read_csv(csv_path, sep="\t")
    # is_AMP can be string "True"/"False" or bool
    n_total = len(df)
    n_amp = (df["is_AMP"] == True) | (df["is_AMP"].astype(str).str.lower() == "true")
    pct = 100 * n_amp.sum() / n_total if n_total else 0
    labels.append(os.path.basename(d).replace("_macrel", ""))
    pct_amp.append(pct)

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(labels, pct_amp, color="steelblue", edgecolor="black")
ax.set_ylabel("% is_AMP")
ax.set_title("Percentage of AMP predictions per dataset")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "figures", "macrel_is_amp_barplot.png"), dpi=150, bbox_inches="tight")
plt.close()