from __future__ import annotations

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "FIG_covalent")

CA2_ROUNDS = [1, 2, 3, 4, 5, 6, 7, 8]
CA2_MARGIN = [0.88, 0.78, 0.75, 0.41, -0.91, -1.07, -2.33, -2.37]
CA2_COVALENT_ROUND = 5

WARHEADS = ["chloroacetamide", "propiolamide", "acrylamide", "acrylonitrile"]
DE_CYS = [-7.25, -38.62, -29.59, -27.99]
DE_SER = [-7.81, -45.51, -29.07, -28.60]

LEADS = [("ABL1 over SRC", "chloroacetamide", "Cys388"),
         ("CA2 over CA1", "propiolamide", "Cys206"),
         ("EGFR T790M", "chloroacetamide", "Cys797"),
         ("KRAS G12C", "acrylonitrile", "Cys12"),
         ("PIK3CA over PIK3CD", "acrylamide", "-")]

def main() -> int:
    S.use()
    fig = plt.figure(figsize=(S.COL_2, 2.35))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.05, 1.15], wspace=0.34,
                          left=0.062, right=0.988, top=0.83, bottom=0.235)

    ax = fig.add_subplot(gs[0])
    ax.axhline(0, color=S.C_RULE, linewidth=0.6, zorder=1)
    ax.plot(CA2_ROUNDS, CA2_MARGIN, marker="o", markersize=3.1, linewidth=1.0,
            color=S.C_GEN, zorder=3)
    i = CA2_ROUNDS.index(CA2_COVALENT_ROUND)
    ax.plot([CA2_COVALENT_ROUND], [CA2_MARGIN[i]], marker="o", markersize=5.4,
            color=S.C_LOSE, zorder=4)
    ax.annotate(S.clean("r5: propiolamide at Cys206\nCA1 has Ser206"),
                xy=(CA2_COVALENT_ROUND, CA2_MARGIN[i]), xycoords="data",
                xytext=(0.34, 0.34), textcoords="axes fraction",
                fontsize=5.4, color=S.C_LOSE, linespacing=1.3, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", linewidth=0.5, color=S.C_LOSE,
                                shrinkA=0, shrinkB=2))
    ax.set_xlabel(S.clean("growth round"))
    ax.set_ylabel(S.clean("selectivity margin (kcal/mol)"))
    ax.set_xticks(CA2_ROUNDS)
    S.despine(ax)
    S.panel_label(ax, "a", dx=-0.26, dy=1.10)
    ax.set_title(S.clean("the score punishes the covalent round"), fontsize=6.6, loc="left",
                 pad=11)

    ax2 = fig.add_subplot(gs[1])
    y = np.arange(len(WARHEADS))
    h = 0.36
    ax2.barh(y + h / 2, DE_CYS, height=h, color=S.C_GEN, label=S.clean("cysteine (Cys-SH)"))
    ax2.barh(y - h / 2, DE_SER, height=h, color=S.C_REF, label=S.clean("serine (Ser-OH)"))
    ax2.set_yticks(y)
    ax2.set_yticklabels([S.clean(w) for w in WARHEADS], fontsize=5.9)
    ax2.set_xlabel(S.clean("reaction energy (kcal/mol)"))
    ax2.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), ncol=2, fontsize=5.4,
               handlelength=1.1, columnspacing=1.2)
    S.despine(ax2)
    S.panel_label(ax2, "b", dx=-0.30, dy=1.10)
    ax2.set_title(S.clean("thermodynamics favours serine"), fontsize=6.6, loc="left", pad=11)
    ax2.text(0.0, 1.015, S.clean("GFN2-xTB, ALPB water, both neutral, same leaving group"),
             transform=ax2.transAxes, fontsize=5.1, ha="left", va="bottom", color="#888888")

    ax3 = fig.add_subplot(gs[2])
    ax3.axis("off")
    S.panel_label(ax3, "c", dx=-0.04, dy=1.10)
    ax3.set_title(S.clean("endpoint energies do not explain it"), fontsize=6.6,
                   loc="left", pad=11)
    txt = ("At pH 7.4 a cysteine side chain (pKa about 8.3, lower in\n"
           "many pockets) is partly THIOLATE, among the strongest\n"
           "nucleophiles in a protein. A serine hydroxyl (pKa about 13)\n"
           "is never ionised and attacks as a neutral alcohol.\n\n"
           "The difference is in the barrier, not the product. Neither a\n"
           "docking score nor an endpoint reaction energy contains a\n"
           "barrier, so covalent rounds are reported separately here\n"
           "and are never folded into the selectivity margin.")
    ax3.text(0.0, 0.94, S.clean(txt), fontsize=5.6, va="top", ha="left", color="#333333",
             linespacing=1.5, transform=ax3.transAxes)
    rows = "\n".join("%-20s %-16s %s" % (c, w, r) for c, w, r in LEADS)
    ax3.text(0.0, 0.245, S.clean("final leads carrying a warhead\n" + rows),
             fontsize=5.2, va="top", ha="left", color=S.C_TEXT, family="monospace",
             linespacing=1.5, transform=ax3.transAxes)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote %s.pdf and .png" % os.path.relpath(OUT, HERE))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
