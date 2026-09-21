from __future__ import annotations

import json
import os
import sys
import textwrap
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S
from render_gui_rounds import describe, divergence_note

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "gui_rounds")
OUT = os.path.join(HERE, "FIG_gui_rounds")

INK, GREY, GREEN, ORANGE = "#111111", "#5a5a5a", "#16834B", "#B4553C"

CROP = (700, 300, 1820, 1080)

def rows() -> List[Dict]:
    man = os.path.join(SRC, "manifest.json")
    if not os.path.exists(man):
        return []
    out, prev = [], None
    for r in json.load(open(man, encoding="utf-8")):
        png = os.path.join(SRC, "view_round%d.png" % r["round"])
        if not os.path.exists(png):
            continue
        r["png"] = png
        r["did"] = describe(prev, r["smiles"])
        r["note"] = divergence_note(r.get("goal_residue") or "")
        out.append(r)
        prev = r["smiles"]
    return out

def main() -> int:
    S.use()
    R = rows()
    if not R:
        print("run build_gui_rounds.py then render_gui_rounds.py first")
        return 1

    n = len(R)
    fig = plt.figure(figsize=(S.COL_2, 3.25))
    gs = fig.add_gridspec(1, n, wspace=0.035, left=0.012, right=0.988,
                          top=0.870, bottom=0.235)

    for i, r in enumerate(R):
        ax = fig.add_subplot(gs[0, i])
        img = mpimg.imread(r["png"])
        l, t, rr, b = CROP
        ax.imshow(img[t:b, l:rr])
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("#d8dcd8")
            s.set_linewidth(0.6)

        last = (i == n - 1)
        ax.set_title(S.clean("%s   round %d   aim %s%s"
                             % ("abc"[i], r["round"], r.get("goal_residue") or "",
                                "   carried lead" if last else "")),
                     loc="left", fontsize=7.2, fontweight="bold",
                     color=GREEN if last else INK, pad=3.5)

        first = ("carried from round %d" % (r["round"] - 1)) if i else "first round"
        ax.text(0.5, -0.045, S.clean(first), transform=ax.transAxes, fontsize=5.8,
                color=GREY, ha="center", va="top")
        ax.text(0.5, -0.115, S.clean(textwrap.fill(r["did"], 38)), transform=ax.transAxes,
                fontsize=5.8, color=INK, ha="center", va="top", linespacing=1.4)
        note = r.get("note") or ""
        ax.text(0.5, -0.255, S.clean(note), transform=ax.transAxes, fontsize=5.8,
                color=(ORANGE if note.startswith("CA1 has") else GREY),
                ha="center", va="top",
                fontweight=("bold" if note.startswith("CA1 has") else "normal"))

    fig.text(0.012, 0.955,
             S.clean("Sequential residue-conditioned design in the binding pocket"),
             fontsize=8.4, fontweight="bold", color=INK)

    for ext in ("pdf", "png", "svg"):
        fig.savefig("%s.%s" % (OUT, ext), dpi=600, facecolor="white",
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf, .png, .svg  (%d panels)" % (os.path.relpath(OUT, HERE), n))
    for r in R:
        print("  round %d  aim %-8s  %-34s  %s"
              % (r["round"], r.get("goal_residue"), r["did"], r.get("note")))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
