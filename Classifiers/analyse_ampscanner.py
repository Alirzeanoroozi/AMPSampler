import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_DIR = "../results"
AMPSCANNER_DIR = os.path.join(RESULTS_DIR, "AMPScanner")

# Find all AMPScanner CSV files
csv_files = sorted(glob.glob(os.path.join(AMPSCANNER_DIR, "*_ampscanner.csv")))

labels = []
pct_amp = []

for csv_path in csv_files:
    df = pd.read_csv(csv_path)
    n_total = len(df)
    n_amp = (df["Prediction_Class"] == "AMP").sum()
    pct = 100 * n_amp / n_total if n_total else 0
    labels.append(os.path.basename(csv_path).replace("_ampscanner.csv", ""))
    pct_amp.append(pct)

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(labels, pct_amp, color="steelblue", edgecolor="black")
ax.set_ylabel("% AMP")
ax.set_title("AMPScanner: Percentage of AMP predictions per dataset")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
os.makedirs(os.path.join(RESULTS_DIR, "figures"), exist_ok=True)
plt.savefig(os.path.join(RESULTS_DIR, "figures", "ampscanner_is_amp_barplot.png"), dpi=150, bbox_inches="tight")
plt.close()