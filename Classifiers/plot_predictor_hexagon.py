#!/usr/bin/env python3
"""Create a polygon/radar plot from the predictor-comparison table."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METRICS = [
    "AMP-Scanner\nAMP Rate",
    "AMP-Scanner\nProb",
    "AMPlify\nAMP Rate",
    "AMPlify\nProb",
    "HydrAMP\nAMP",
    "HydrAMP\nMIC",
    "MACREL\nAMP",
    "MACREL\nHemo",
]

DATA = {
    "AMPGAN": [0.6325, 0.6311, 0.6541, 0.6350, 0.6314, 0.3163, 0.4744, 0.5978],
    "PepCVAE": [0.4954, 0.5046, 0.4917, 0.4913, 0.5155, 0.2052, 0.3612, 0.4106],
    "HydrAMP": [0.6392, 0.6373, 0.6053, 0.5859, 0.7818, 0.5120, 0.3976, 0.4726],
    "AMP-Diffusion": [0.6948, 0.6953, 0.1242, 0.1664, 0.8147, 0.4560, 0.4465, 0.5001],
    "DiT-AMP (Ours)": [0.9980, 0.9947, 0.9970, 0.9885, 0.9748, 0.7953, 0.8110, 0.9423],
}


def normalize(values: np.ndarray) -> np.ndarray:
    """Min-max normalize each metric, inverting lower-is-better metrics."""
    normalized = np.zeros_like(values, dtype=float)

    for idx, metric in enumerate(METRICS):
        column = values[:, idx]
        min_value = column.min()
        max_value = column.max()

        if np.isclose(max_value, min_value):
            normalized[:, idx] = 1.0
            continue

        scaled = (column - min_value) / (max_value - min_value)
        normalized[:, idx] = scaled

    return normalized


def main() -> None:
    labels = list(DATA)
    values = np.array([DATA[label] for label in labels], dtype=float)
    # values = normalize(raw_values)

    angles = np.linspace(0, 2 * np.pi, len(METRICS), endpoint=False)
    closed_angles = np.concatenate([angles, angles[:1]])

    fig, ax = plt.subplots(figsize=(8.5, 8.5), subplot_kw={"projection": "polar"})

    for label, row in zip(labels, values):
        closed_row = np.concatenate([row, row[:1]])
        linewidth = 2.6 if label == "DiT-AMP (Ours)" else 1.6
        alpha = 0.18 if label == "DiT-AMP (Ours)" else 0.08
        ax.plot(closed_angles, closed_row, linewidth=linewidth, label=label)
        ax.fill(closed_angles, closed_row, alpha=alpha)

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(METRICS, fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=8)
    ax.grid(alpha=0.35)
    ax.set_title("Normalized Predictor Metrics", pad=24, fontsize=14)
    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.12), fontsize=9)

    output_path = Path(__file__).resolve().parent / "images" / "predictor_hexagon_plot.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
