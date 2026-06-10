import os
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../data/tables/boltz2_confidence_summary_NDM5_KPC3_boltz2.csv")

fig, ax = plt.subplots(figsize=(8, 6))
for target in df["target"].unique():
    subset = df[df["target"] == target]
    ax.scatter(
        subset["confidence_score"],
        subset["iptm"],
        label=target,
        alpha=0.7,
    )
ax.set_xlabel("Confidence")
ax.set_ylabel("iptm")
ax.set_title("iptm vs Confidence by target")
ax.legend()
plt.tight_layout()
os.makedirs("../results/figures", exist_ok=True)
plt.savefig("../results/figures/boltz2_iptm_confidence.png", dpi=150, bbox_inches="tight")
plt.close()
