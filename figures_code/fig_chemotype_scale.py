from __future__ import annotations

import csv
import os
import sys
from typing import Dict, List

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
SCALE_CSV = os.path.join(os.path.dirname(HERE), "results", "STATS_chemotype_scale.csv")
OUT = os.path.join(HERE, "FIG_chemotype_scale")

C_LEAD = "#1b6ca8"
C_WIN = "#2a7d4f"
C_NULL = "#b0b7bd"
C_SCALE = "#8c8c8c"
C_WARN = "#b4472e"

def pair_scale() -> Dict[str, List[float]]:
    import itertools
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from stats_structural import mcs_pct, load

    _, drugs, _ = load()
    out: Dict[str, List[float]] = {}
    for case, ds in drugs.items():
        mols = [(d["key"], Chem.MolFromSmiles(d["smiles"])) for d in ds]
        mols = [(n, m) for n, m in mols if m is not None]
        if len(mols) < 2:
            continue
        out[case] = sorted(mcs_pct(a, b) for (_, a), (_, b) in itertools.combinations(mols, 2))
    return out

def rows() -> List[Dict]:
    if not os.path.exists(SCALE_CSV):
        return []
    out = []
    for r in csv.DictReader(open(SCALE_CSV, newline="", encoding="utf-8")):
        try:
            r["lead_best"] = float(r["lead_best"])
            r["null_median"] = float(r["null_median"])
            r["drug_to_drug_median"] = (float(r["drug_to_drug_median"])
                                        if r.get("drug_to_drug_median") else None)
            r["p"] = float(r["p"])
        except (TypeError, ValueError):
            continue
        out.append(r)
    out.sort(key=lambda r: -r["lead_best"])
    return out

def main() -> int:
    S.use()
    R = rows()
    if not R:
        print("run stats_chemotype_scale.py first")
        return 1
    scale = pair_scale()
    allpairs = sorted(v for vs in scale.values() for v in vs)
    med_all = float(np.median(allpairs)) if allpairs else float("nan")

    fig = plt.figure(figsize=(S.COL_2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.22], wspace=0.30,
                          left=0.085, right=0.985, top=0.90, bottom=0.175)

    axl = fig.add_subplot(gs[0, 0])
    seen_drugsets = set()
    uniq = []
    for c in scale:
        key = tuple(round(v, 3) for v in scale[c])
        if key in seen_drugsets:
            continue
        seen_drugsets.add(key)
        uniq.append(c)
    order = sorted(uniq, key=lambda c: -float(np.median(scale[c])))
    ys = list(range(len(order)))[::-1]
    for y, case in zip(ys, order):
        v = scale[case]
        axl.plot([min(v), max(v)], [y, y], linewidth=0.8, color=C_NULL, zorder=1)
        axl.scatter(v, [y] * len(v), s=6, color=C_SCALE, zorder=2, linewidths=0)
        axl.scatter([float(np.median(v))], [y], s=17, color="#333333", zorder=3,
                    marker="|", linewidths=1.2)
    axl.axvline(med_all, color=C_SCALE, linewidth=0.8, linestyle=(0, (4, 2)), zorder=0)
    axl.set_yticks(ys)
    axl.set_yticklabels([c.split("_")[0][:12] for c in order], fontsize=5.2)
    axl.set_xlabel(S.clean("shared core between two marketed drugs, %"), fontsize=5.9)
    axl.tick_params(axis="x", labelsize=5.2)
    axl.set_xlim(0, 100)
    for sp in ("top", "right"):
        axl.spines[sp].set_visible(False)
    axl.set_title(S.clean("What two real drugs for one target share"),
                  fontsize=7.0, loc="left", color=S.C_TEXT, pad=5)

    axr = fig.add_subplot(gs[0, 1])
    ys = list(range(len(R)))[::-1]
    for y, r in zip(ys, R):
        beats = (r["drug_to_drug_median"] is not None
                 and r["lead_best"] > r["drug_to_drug_median"])
        col = C_WIN if beats else C_LEAD
        axr.plot([r["null_median"], r["lead_best"]], [y, y], linewidth=0.7,
                 color=("#cfd6dc" if not beats else "#bcd6c6"), zorder=1)
        axr.scatter([r["null_median"]], [y], s=13, color=C_NULL, zorder=2, linewidths=0)
        if r["drug_to_drug_median"] is not None:
            axr.scatter([r["drug_to_drug_median"]], [y], s=26, facecolors="none",
                        edgecolors=C_SCALE, linewidths=0.8, zorder=3)
        axr.scatter([r["lead_best"]], [y], s=26, color=col, zorder=4, linewidths=0)
        lab = "%.0f%%" % r["lead_best"]
        if r["p"] < 0.05:
            lab += "   p=%.3f" % r["p"]
        axr.text(99.0, y, S.clean(lab), fontsize=4.9, va="center", ha="right",
                 color=(col if r["p"] < 0.05 else "#555555"),
                 fontweight=("bold" if r["p"] < 0.05 else "normal"))

    axr.set_yticks(ys)
    axr.set_yticklabels(["%s  (%s)" % (r["case"].split("_")[0][:9],
                                       (r["matched_drug"] or "")[:11]) for r in R], fontsize=5.2)
    axr.set_xlabel(S.clean("best shared core with any marketed drug for that target, %"),
                   fontsize=5.9)
    axr.tick_params(axis="x", labelsize=5.2)
    axr.set_xlim(0, 100)
    for sp in ("top", "right"):
        axr.spines[sp].set_visible(False)
    axr.set_title(S.clean("Where the generated leads fall on it"),
                  fontsize=7.0, loc="left", color=S.C_TEXT, pad=5)

    n_beat = sum(1 for r in R if r["drug_to_drug_median"] is not None
                 and r["lead_best"] > r["drug_to_drug_median"])
    n_sig = sum(1 for r in R if r["p"] < 0.05)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), dpi=500 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf and .png" % os.path.relpath(OUT, HERE))
    print("  drug-to-drug median across %d pairs: %.0f%%" % (len(allpairs), med_all))
    print("  %d of %d leads above their target's drug-to-drug median" % (n_beat, len(R)))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
