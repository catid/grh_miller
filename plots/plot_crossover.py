#!/usr/bin/env python3
"""Render the tweet-sized latency/status diagram from crossover_results.csv."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "crossover_results.csv"
OUT = ROOT / "tweet_performance_diagram.png"

COLORS = {
    "aks_textbook": "#D1495B",             # red: exact
    "miller_grh_reference": "#00798C",     # teal: conditional
    "gmpy2_probable_prime": "#EDAE49",     # gold: probable
}
LABELS = {
    "aks_textbook": "AKS textbook  |  exact / unconditional",
    "miller_grh_reference": "Miller GRH scan  |  conditional (ERH)",
    "gmpy2_probable_prime": "GMP probable-prime  |  probabilistic filter",
}
ORDER = ["aks_textbook", "miller_grh_reference", "gmpy2_probable_prime"]


def read_data() -> dict[str, list[dict[str, str]]]:
    # Allow provenance comments in the compact CSV while keeping it normal CSV.
    with DATA.open(newline="") as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    out = {name: [] for name in ORDER}
    for row in rows:
        row["bits"] = int(row["bits"])
        row["timeout_seconds"] = float(row["timeout_seconds"])
        row["median_seconds"] = float(row["median_seconds"]) if row["median_seconds"] else None
        out[row["algorithm"]].append(row)
    return out


def plot() -> Path:
    data = read_data()
    fig, ax = plt.subplots(figsize=(12, 6.75), dpi=100)
    fig.patch.set_facecolor("#FBFCFE")
    ax.set_facecolor("#FBFCFE")
    for name in ORDER:
        rows = sorted(data[name], key=lambda r: r["bits"])
        complete = [r for r in rows if r["observed"] == "completed"]
        censored = [r for r in rows if r["observed"] == "timeout"]
        color = COLORS[name]
        ax.plot(
            [r["bits"] for r in complete],
            [r["median_seconds"] for r in complete],
            marker="o", markersize=5, linewidth=2.5, color=color,
            label=LABELS[name], zorder=3,
        )
        # Censored measurements are plotted exactly at the threshold with an x;
        # they mean “> 0.5 s”, never an observed 0.5-second latency.
        if censored:
            threshold = censored[0]["timeout_seconds"]
            ax.scatter(
                [r["bits"] for r in censored], [threshold] * len(censored),
                marker="x", s=56, linewidths=2.0, color=color, zorder=4,
            )

    timeout = 0.5
    ax.axhline(timeout, color="#5B6573", linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
    ax.text(8.6, timeout * 1.12, "timeout ceiling 0.5 s  (× = > ceiling)",
            color="#4A5563", fontsize=8.5, va="bottom", ha="left")
    ax.annotate(
        "AKS times out\nfrom 20 bits",
        xy=(20, timeout), xytext=(24, 0.16), color=COLORS["aks_textbook"],
        fontsize=8.5, fontweight="bold", ha="left", va="top",
        arrowprops=dict(arrowstyle="->", color=COLORS["aks_textbook"], lw=1.2),
    )
    ax.annotate(
        "GRH scan times out\nfrom 256 bits",
        xy=(256, timeout), xytext=(160, 0.16), color=COLORS["miller_grh_reference"],
        fontsize=8.5, fontweight="bold", ha="left", va="top",
        arrowprops=dict(arrowstyle="->", color=COLORS["miller_grh_reference"], lw=1.2),
    )

    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    xticks = [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(x) for x in xticks])
    ax.set_xlim(7, 9500)
    ax.set_ylim(1e-7, 1.15)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=9))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", axis="both", color="#B8C1CC", alpha=0.34, linewidth=0.8)
    ax.grid(which="minor", axis="y", color="#D7DDE5", alpha=0.22, linewidth=0.5)
    ax.set_xlabel("input size (bits)", fontsize=11)
    ax.set_ylabel("median latency per input (seconds, log scale)", fontsize=11)
    ax.set_title("")
    fig.text(
        0.085, 0.955, "Exact vs conditional vs probabilistic primality checks",
        ha="left", va="top", fontsize=18, fontweight="bold", color="#17202A",
    )
    fig.text(
        0.085, 0.918,
        "Constructed-prime corpus  |  lower is faster  |  censored timeouts are not 0.5 s measurements",
        ha="left", va="top", fontsize=10, color="#4A5563",
    )
    leg = ax.legend(
        loc="lower right", fontsize=9.5, frameon=True, framealpha=0.97,
        facecolor="white", edgecolor="#C8D0DA", borderpad=0.65,
    )
    for line in leg.get_lines():
        line.set_linewidth(3)

    # Explicit status key makes the proof boundary readable in a screenshot.
    cards = [
        ("EXACT", "AKS: unconditional algorithm", COLORS["aks_textbook"]),
        ("CONDITIONAL", "Miller: ERH / B(n)=ceil(2 ln(n)^2)", COLORS["miller_grh_reference"]),
        ("PROBABLE", "GMP: probable-prime filter; no proof", COLORS["gmpy2_probable_prime"]),
    ]
    x0, y0, width, height, gap = 0.02, -0.235, 0.30, 0.095, 0.02
    for i, (head, body, color) in enumerate(cards):
        x = x0 + i * (width + gap)
        ax.text(x, y0 + height * 0.64, head, transform=ax.transAxes,
                fontsize=9, fontweight="bold", color="white", va="center",
                bbox=dict(boxstyle="round,pad=0.32", facecolor=color, edgecolor=color))
        ax.text(x + 0.105, y0 + height * 0.64, body, transform=ax.transAxes,
                fontsize=8.5, color="#27313B", va="center")

    fig.text(
        0.02, 0.018,
        "Source: results/crossover_results.csv (derived from the frozen targeted run, 2026-11-11) | "
        "CPU 5, 3 repeats (2 at 4096+), 0.5 s timeout | Python 3.13.13 + gmpy2 2.3.2/GMP 6.3.0",
        fontsize=7.5, color="#5B6573", ha="left",
    )
    fig.subplots_adjust(left=0.085, right=0.985, top=0.82, bottom=0.31)
    fig.savefig(OUT, dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)
    return OUT


if __name__ == "__main__":
    path = plot()
    print(f"wrote {path} ({path.stat().st_size} bytes)")
