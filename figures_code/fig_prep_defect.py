from __future__ import annotations

import os
import sys
import urllib.request
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_pdb_cache")
OUT = os.path.join(HERE, "FIG_prep_defect")

C_BAD = "#b4472e"
C_OK = "#2a7d4f"
C_GREY = "#8c8c8c"

RECEPTORS = ["6SBL", "7Q0D", "6NPV", "7NG7", "6TFV", "8A27", "6OIM", "9IAY",
             "7NH5", "9C1W", "7AAC", "7R59"]

def fetch(pdb_id: str) -> Optional[str]:
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "%s.pdb" % pdb_id.upper())
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return open(path, encoding="utf-8", errors="ignore").read()
    try:
        txt = urllib.request.urlopen(
            "https://files.rcsb.org/download/%s.pdb" % pdb_id.upper(), timeout=60
        ).read().decode("utf-8", "ignore")
    except Exception:
        return None
    open(path, "w", encoding="utf-8").write(txt)
    return txt

def damage(pdb_text: str):
    corrupted, altloc, example = 0, 0, []
    residues = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or len(line) < 20:
            continue
        if line[16] == " ":
            continue
        altloc += 1
        broken = line.replace(line[16], " ", 1)
        fixed = line[:16] + " " + line[17:]
        if broken != fixed:
            corrupted += 1
            residues.add((line[21], line[22:26].strip()))
            if not example:
                example = [line, broken, fixed]
    return corrupted, altloc, example, len(residues)

def main() -> int:
    S.use()
    rows = []
    example = []
    for pid in RECEPTORS:
        txt = fetch(pid)
        if not txt:
            continue
        c, a, ex, nres = damage(txt)
        rows.append({"pdb": pid, "corrupted": c, "altloc": a, "residues": nres})
        if ex and (not example or pid == "6SBL"):
            example = ex
    if not rows:
        print("no receptors could be retrieved")
        return 1
    rows.sort(key=lambda r: -r["corrupted"])

    fig = plt.figure(figsize=(S.COL_2, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.32, 1.0], wspace=0.18,
                          left=0.015, right=0.985, top=0.90, bottom=0.10)

    axm = fig.add_subplot(gs[0, 0])
    axm.axis("off"); axm.set_xlim(0, 1); axm.set_ylim(0, 1)
    axm.text(0.0, 1.0, S.clean("The mechanism"), fontsize=7.4, fontweight="bold",
             color=S.C_TEXT, va="top")
    axm.text(0.0, 0.925,
             S.clean("str.replace takes a character, not a column. The atom name holds the "
                     "altloc letter first."),
             fontsize=5.3, color=C_GREY, va="top")

    if example:
        src, broken, fixed = example[0], example[1], example[2]
        seg = slice(0, 30)
        axm.text(0.0, 0.80, S.clean("deposited record"), fontsize=5.2,
                 color=S.C_TEXT, va="top", fontweight="bold")
        axm.text(0.0, 0.735, src[seg], fontsize=5.6, family="monospace", va="top",
                 color=S.C_TEXT)
        axm.text(0.0, 0.665, " " * 13 + "^" * 4 + " " + "^", fontsize=5.6, family="monospace",
                 va="top", color=C_BAD)
        axm.text(0.0, 0.615, S.clean("             atom name    altloc, column 17"),
                 fontsize=4.7, family="monospace", va="top", color=C_BAD)

        axm.text(0.0, 0.50, S.clean("after the broken transform"), fontsize=5.2,
                 color=C_BAD, va="top", fontweight="bold")
        axm.text(0.0, 0.435, broken[seg], fontsize=5.6, family="monospace", va="top",
                 color=C_BAD)
        axm.text(0.0, 0.365,
                 S.clean("the letter was removed from the ATOM NAME; the altloc is still there"),
                 fontsize=4.7, color=C_BAD, va="top")

        axm.text(0.0, 0.25, S.clean("after the corrected transform"), fontsize=5.2,
                 color=C_OK, va="top", fontweight="bold")
        axm.text(0.0, 0.185, fixed[seg], fontsize=5.6, family="monospace", va="top", color=C_OK)
        axm.text(0.0, 0.115, S.clean("the atom keeps its name; column 17 is blanked"),
                 fontsize=4.7, color=C_OK, va="top")

    axm.text(0.0, 0.02,
             S.clean("Consequence: the residue carries two atoms with the same name about 2.5 A "
                     "apart, so the pocket gains contacts it does not have."),
             fontsize=5.0, color=S.C_TEXT, va="top")

    axb = fig.add_subplot(gs[0, 1])
    ys = list(range(len(rows)))[::-1]
    vals = [r["corrupted"] for r in rows]
    axb.barh(ys, vals, height=0.62,
             color=[C_BAD if v else C_GREY for v in vals], edgecolor="none")
    axb.set_yticks(ys)
    axb.set_yticklabels([r["pdb"] for r in rows], fontsize=5.6)
    axb.set_xlabel(S.clean("ATOM records corrupted"), fontsize=5.8)
    axb.tick_params(axis="x", labelsize=5.4)
    for sp in ("top", "right"):
        axb.spines[sp].set_visible(False)
    for y, v in zip(ys, vals):
        axb.text(v + max(vals) * 0.015, y, S.clean(str(v)), fontsize=5.2, va="center",
                 color=(C_BAD if v else C_GREY))
    axb.set_title(S.clean("Every receptor in the panel"), fontsize=7.0, loc="left",
                  color=S.C_TEXT, pad=4)

    n_hit = sum(1 for r in rows if r["corrupted"])

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), dpi=500 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf and .png" % os.path.relpath(OUT, HERE))
    for r in rows:
        print("  %-6s %4d corrupted of %4d altloc records, %3d residues"
              % (r["pdb"], r["corrupted"], r["altloc"], r["residues"]))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
