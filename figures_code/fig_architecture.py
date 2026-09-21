from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "FIG_architecture")

CASE = "CA2 over CA1"
ROUNDS = [
    ("scaffold", "36 known actives", "benzenesulfonamide", False),
    ("r1", "His94", "zinc anchor", False),
    ("r2", "Phe131", "aromatic", False),
    ("r3", "Leu204", "hydrophobic", False),
    ("r4", "Thr200", "H-bond donor", False),
    ("r5", "Cys206", "covalent", True),
]

def main() -> int:
    S.use()
    fig = plt.figure(figsize=(S.COL_2, 3.35))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.3, 1.45], hspace=0.55,
                          left=0.045, right=0.985, top=0.93, bottom=0.07)

    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, 10); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0, 0.92, S.clean("a   Containment: each round must contain the last"),
            fontsize=7.5, fontweight="bold", va="top", color=S.C_TEXT)
    n = len(ROUNDS)
    for i, (tag, res, feat, cov) in enumerate(ROUNDS):
        x = 0.15 + i * (9.6 / n)
        w = 9.6 / n - 0.35
        ax.add_patch(FancyBboxPatch((x, 0.20), w, 0.36,
                                    boxstyle="round,pad=0.012,rounding_size=0.03",
                                    linewidth=0.7,
                                    edgecolor=(S.C_WIN if cov else "#555555"),
                                    facecolor=("#eef6f0" if cov else "#f6f7f8")))
        block = "%s\n%s" % (res, feat)
        ax.text(x + w / 2, 0.38, S.clean(block), ha="center", va="center",
                fontsize=5.4, color="#444444", linespacing=1.6)
        ax.text(x + w / 2, 0.60, S.clean(tag), ha="center", va="bottom", fontsize=6.2,
                fontweight="bold", color=(S.C_WIN if cov else S.C_TEXT))
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.02, 0.38), (x + w + 0.30, 0.38),
                                         arrowstyle="-|>", mutation_scale=6,
                                         linewidth=0.6, color="#888888"))
    ax.text(0.15, 0.045, S.clean("every box contains the one before it, checked as a subgraph "
                                "match before any candidate is scored"),
            fontsize=5.7, color="#777777", va="center")

    ax2 = fig.add_subplot(gs[1])
    ax2.set_xlim(0, 10); ax2.set_ylim(0, 1); ax2.axis("off")
    ax2.text(0, 0.96, S.clean("b   Residue conditioning: the target of each round is chosen, "
                              "not sampled"),
             fontsize=7.5, fontweight="bold", va="top", color=S.C_TEXT)
    steps = [
        ("target, anti-target", "two structures of\nthe same fold"),
        ("sequence alignment", "which pocket positions\nDIFFER"),
        ("reach test", "which of those a single\ngroup can actually touch"),
        ("model chooses", "one residue, one feature\ntype, with its reasoning"),
        ("engine enumerates", "the chemistry that\nsatisfies it"),
        ("docked and gated", "containment, pose,\nPareto front"),
    ]
    for i, (head, body) in enumerate(steps):
        x = 0.15 + i * 1.62
        ax2.add_patch(FancyBboxPatch((x, 0.30), 1.42, 0.46,
                                     boxstyle="round,pad=0.012,rounding_size=0.03",
                                     linewidth=0.7, edgecolor="#555555", facecolor="white"))
        ax2.text(x + 0.71, 0.68, S.clean(head), ha="center", fontsize=6.0,
                 fontweight="bold", color=S.C_TEXT)
        ax2.text(x + 0.71, 0.50, S.clean(body), ha="center", va="center", fontsize=5.5,
                 color="#555555", linespacing=1.35)
        if i < len(steps) - 1:
            ax2.add_patch(FancyArrowPatch((x + 1.44, 0.53), (x + 1.60, 0.53),
                                          arrowstyle="-|>", mutation_scale=6,
                                          linewidth=0.6, color="#888888"))
    ax2.text(0.15, 0.17, S.clean("the alignment supplies the hypothesis; the model states which "
                                 "one to test this round and why; the engine and the docking "
                                 "decide whether it worked"),
             fontsize=5.7, color="#777777", va="center")

    ax3 = fig.add_subplot(gs[2])
    ax3.set_xlim(0, 10); ax3.set_ylim(0, 1); ax3.axis("off")
    ax3.text(0, 1.12, S.clean("c   What the scoring function measures, and what it does not"),
             fontsize=7.5, fontweight="bold", va="top", color=S.C_TEXT)
    left = ("MEASURED by docking\n"
            "non-covalent contacts in the target pocket\n"
            "the same pose scored in a chimeric counter-pocket\n"
            "the difference between them, the selectivity margin")
    right = ("NOT MEASURED by docking\n"
             "the covalent bond a warhead forms\n"
             "that a thiolate (Cys, pKa 8.3) is a far better nucleophile\n"
             "than an alcohol (Ser, pKa 13) at pH 7.4\n"
             "reaction energies do NOT supply the missing term: balanced\n"
             "like for like they favour serine by up to 7 kcal/mol, so the\n"
             "selectivity is kinetic and covalent rounds are reported apart")
    ax3.add_patch(FancyBboxPatch((0.15, -0.06), 4.6, 0.92,
                                 boxstyle="round,pad=0.014,rounding_size=0.03",
                                 linewidth=0.7, edgecolor="#555555", facecolor="#f6f7f8"))
    ax3.add_patch(FancyBboxPatch((5.05, -0.06), 4.75, 0.92,
                                 boxstyle="round,pad=0.014,rounding_size=0.03",
                                 linewidth=0.7, edgecolor=S.C_WIN, facecolor="#eef6f0"))
    ax3.text(0.30, 0.82, S.clean(left), fontsize=4.9, va="top", color="#444444",
             linespacing=1.5)
    ax3.text(5.20, 0.82, S.clean(right), fontsize=4.9, va="top", color="#2a5c3a",
             linespacing=1.5)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote %s.pdf and .png" % os.path.relpath(OUT, HERE))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
