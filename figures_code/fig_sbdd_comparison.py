import csv
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
from rdkit import Chem, RDLogger, RDConfig
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import QED, Descriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

plt.style.use(["science", "nature", "no-latex"])
HERE = os.path.dirname(os.path.abspath(__file__))
LIVE = ("panel_plan_", "panel_build_")

INK = "#111111"
ACC = "#A63603"
PUB = "#6e6e6e"
REF = "#08519C"

PUBLISHED = [
    ("Reference ligand", -6.36, -7.45, 0.48, 0.73, 22.8, "Luo T1 / TargetDiff T1 / MolCRAFT T2"),
    ("liGAN",            None,  -6.33, 0.39, 0.59, None, "TargetDiff T1 re-eval"),
    ("AR (3D-SBDD)",    -6.344, None,  0.525, 0.657, None, "Luo T1"),
    ("Pocket2Mol (QVina)", -7.288, None, 0.563, 0.765, 17.7, "Pocket2Mol T1 (QVina)"),
    ("TargetDiff",      -5.47,  -7.80, 0.48, 0.58, 24.2, "TargetDiff T1"),
    ("DecompDiff",      -5.67,  -8.39, 0.45, 0.61, 29.4, "DecompDiff T3"),
    ("MolCRAFT",        -6.59,  -7.92, 0.50, 0.69, 22.7, "MolCRAFT T2"),
]
VALIDITY3D = [("DiffSBDD", 0.4), ("liGAN", 2.0), ("Pocket2Mol", 9.0),
              ("3D-SBDD", 10.0), ("TargetDiff", 11.0), ("ResGen", 11.0)]
STRAIN_PUB = [("CrossDocked test", 102.5), ("Pocket2Mol", 194.9), ("3D-SBDD", 592.2),
              ("TargetDiff", 1241.7), ("DiffSBDD", 1243.1)]
RANDOM_FLOOR = 30.83

_pf = FilterCatalogParams(); _pf.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
PAINS = FilterCatalog.FilterCatalog(_pf)
REAGENT = Chem.MolFromSmarts("[N,O,S]-[Cl,Br,I]")
AMIDE = Chem.MolFromSmarts("C(=O)N")
MICHAEL = Chem.MolFromSmarts("C=CC(=O)[N,O]")

def passes_screen(m):
    if m.HasSubstructMatch(REAGENT):
        return False
    if len(m.GetSubstructMatches(MICHAEL)) >= 2:
        return False
    if len(m.GetSubstructMatches(AMIDE)) >= 4:
        return False
    if Descriptors.MolWt(m) > 900:
        return False
    return not PAINS.HasMatch(m)

def octara():
    rows, seen = [], set()
    with open(os.path.join(HERE, "master_results.csv"), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if (r.get("kind") or "") != "generated":
                continue
            if not any((r.get("source_file") or "").startswith(p) for p in LIVE):
                continue
            smi = (r.get("smiles") or "").strip()
            m = Chem.MolFromSmiles(smi) if smi else None
            if m is None or smi in seen or not passes_screen(m):
                continue

            def f(k):
                try:
                    return float((r.get(k) or "").strip())
                except (TypeError, ValueError):
                    return None
            dg = f("affinity_target")
            if dg is None:
                continue
            seen.add(smi)
            rows.append({"dg": dg, "ha": m.GetNumHeavyAtoms(), "qed": QED.qed(m),
                         "sa_raw": sascorer.calculateScore(m), "strain": f("strain_kcal"),
                         "pv": (r.get("pose_valid") or "").strip().lower()})
    return rows

def main():
    R = octara()
    sa_norm = [(10.0 - v["sa_raw"]) / 9.0 for v in R]
    strain = [v["strain"] for v in R if v["strain"] is not None]
    pv = [v["pv"] for v in R if v["pv"] in ("true", "false")]

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.34,
                          left=0.115, right=0.975, top=0.945, bottom=0.115)

    ax = fig.add_subplot(gs[0, :])
    names = [p[0] for p in PUBLISHED if p[1] is not None]
    vals = [p[1] for p in PUBLISHED if p[1] is not None]
    small = [v for v in R if v["ha"] <= 25]
    names.append("Octara ($\\leq$25 HA)")
    vals.append(statistics.mean([v["dg"] for v in small]))
    names.append("Octara (all)")
    vals.append(statistics.mean([v["dg"] for v in R]))
    cols = [REF if n.startswith("Reference") else (ACC if "Octara" in n else PUB) for n in names]
    y = np.arange(len(names))
    ax.barh(y, vals, color=cols, height=0.62, edgecolor=INK, linewidth=0.4)
    ax.set_xlim(min(vals) - 1.35, 0.35)
    for i, v in enumerate(vals):
        ax.text(v - 0.14, i, "%.2f" % v, va="center", ha="right", fontsize=6.4, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.8)
    ax.invert_yaxis()
    ax.set_xlabel("in-place Vina score (kcal mol$^{-1}$), lower is better", fontsize=7.2)
    ax.set_title("A   In-place Vina score", loc="left", fontsize=8.2, fontweight="bold")
    ax.axvline(-6.36, color=REF, lw=0.7, ls=(0, (3, 2)))
    ax.annotate("CrossDocked reference", xy=(-6.36, -0.72), xytext=(-6.36, -0.72),
                fontsize=6, color=REF, ha="center", va="bottom", annotation_clip=False)

    ax = fig.add_subplot(gs[1, 0])
    bins = [(0, 25, "$\\leq$25"), (26, 30, "26-30"), (31, 99, "$>$30")]
    xs, ys, ns = [], [], []
    for lo, hi, lab in bins:
        sub = [v["dg"] for v in R if lo <= v["ha"] <= hi]
        if sub:
            xs.append(lab)
            ys.append(statistics.mean(sub))
            ns.append(len(sub))
    ax.bar(xs, ys, color=ACC, edgecolor=INK, linewidth=0.4, width=0.6)
    ax.set_ylim(min(ys) - 1.6, 0)
    for i, (v, n) in enumerate(zip(ys, ns)):
        ax.text(i, v + 0.30, "n=%d" % n, ha="center", va="bottom", fontsize=6.2, color="white")
    ax.axhline(-6.36, color=REF, lw=0.8, ls=(0, (3, 2)))
    ax.text(-0.42, -6.36, "reference", fontsize=6, color=REF, ha="left", va="bottom")
    ax.set_ylabel("mean Vina (kcal mol$^{-1}$)", fontsize=7.2)
    ax.set_xlabel("heavy atoms", fontsize=7.2)
    ax.set_title("B   Octara Vina score\n     stratified by size", loc="left", fontsize=8.2,
                 fontweight="bold")

    ax = fig.add_subplot(gs[1, 1])
    pairs = sorted([(p[0], p[5], p[4]) for p in PUBLISHED
                    if p[4] is not None and p[5] is not None], key=lambda t: t[1])
    names, theirs, ours, notes = [], [], [], []
    for nm, pha, psa in pairs:
        win = [s for s, v in zip(sa_norm, R) if abs(v["ha"] - pha) <= 3.0]
        names.append("%s\n%.1f HA" % (nm, pha))
        theirs.append(psa)
        ours.append(statistics.mean(win) if win else None)
        notes.append(len(win))
    y = np.arange(len(names))
    h = 0.36
    ax.barh(y + h / 2, theirs, height=h, color=PUB, edgecolor=INK, linewidth=0.4,
            label="published, as reported")
    ax.barh([yy - h / 2 for yy, o in zip(y, ours) if o is not None],
            [o for o in ours if o is not None], height=h, color=ACC,
            edgecolor=INK, linewidth=0.4, label="this study, same size")
    for yy, o, n, tv in zip(y, ours, notes, theirs):
        if o is None:
            ax.text(0.508, yy - h / 2, "no molecules of this size in this study",
                    fontsize=6, color="#8a4a4a", va="center")
        else:
            ax.text(o + 0.006, yy - h / 2, "%.3f  (n=%d)" % (o, n), fontsize=6,
                    color=ACC, va="center")
        ax.text(tv + 0.006, yy + h / 2, "%.3f" % tv, fontsize=6, color="#4f4f4f",
                va="center")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_ylim(len(names) - 0.45, -0.70)
    ax.set_xlim(0.5, 0.97)
    ax.set_xlabel("SA, normalised $(10-\\mathrm{SA})/9$, higher is better", fontsize=6.8)
    ax.set_title("C   Synthetic accessibility, matched by size", loc="left",
                 fontsize=8.2, fontweight="bold")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
              frameon=False, handlelength=1.2, columnspacing=1.4)

    ax = fig.add_subplot(gs[2, 0])
    octara_pv = 100.0 * pv.count("true") / len(pv) if pv else 0.0
    _pd_names = ["MolPilot\n2025", "Octara"]
    _pd_vals = [95.9, octara_pv]
    _pd_cols = ["#8c8c8c", ACC]
    ax.bar(range(len(_pd_vals)), _pd_vals, color=_pd_cols, edgecolor=INK,
           linewidth=0.4, width=0.46)
    for _i, _v in enumerate(_pd_vals):
        ax.text(_i, _v + 1.6, "%.1f%%" % _v, ha="center", fontsize=6.4, color=INK)
    ax.set_xlim(-0.65, len(_pd_vals) - 0.35)
    ax.set_ylim(0, 108)
    ax.set_xticks(range(len(_pd_vals)))
    ax.set_xticklabels(_pd_names, fontsize=6.4)
    ax.set_ylabel("PoseBusters pass rate (\\%)", fontsize=7.2)
    ax.set_title("D   Pose validity, PoseBusters", loc="left", fontsize=8.2,
                 fontweight="bold")

    ax = fig.add_subplot(gs[2, 1])
    med = statistics.median(strain) if strain else 0.0
    ax.bar([0], [med], color=ACC, edgecolor=INK, linewidth=0.4, width=0.42)
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(0, max(1.0, med * 1.45))
    ax.set_xticks([0])
    ax.set_xticklabels(["Octara"], fontsize=6.4)
    ax.text(0, med * 1.06, "%.1f" % med, ha="center", fontsize=6.4, color=INK)
    ax.set_ylabel("median strain (kcal mol$^{-1}$)", fontsize=7.2)
    ax.set_title("E   Conformer strain (this study)", loc="left", fontsize=8.2,
                 fontweight="bold")

    HERE_DIR = os.path.dirname(os.path.abspath(__file__))
    for OUT, stem in ((r"C:\Users\amine\Downloads\octara_figures",
                       "Figure4_vs_published_generators"),
                      (HERE_DIR, "FIG_sbdd_comparison")):
        os.makedirs(OUT, exist_ok=True)
        for ext in ("pdf", "png", "svg"):
            p = os.path.join(OUT, "%s.%s" % (stem, ext))
            if os.path.exists(p):
                os.unlink(p)
            fig.savefig(p, dpi=600, facecolor="white")
    plt.close(fig)

    print("Octara, tonight's engine: %d molecules" % len(R))
    print("  mean Vina all        %.2f" % statistics.mean([v["dg"] for v in R]))
    print("  mean Vina <=25 HA    %.2f  (n=%d)"
          % (statistics.mean([v["dg"] for v in small]), len(small)))
    print("  SA normalised        %.3f  (best published 0.765)" % statistics.mean(sa_norm))
    print("  median strain        %.1f kcal/mol, this study"
          % (statistics.median(strain) if strain else float("nan")))
    print("    published context, NOT ranked against the above (GenBench3D MMFF94s medians): %s"
          % ", ".join("%s %.1f" % (n, v) for n, v in STRAIN_PUB))
    print("  pose validity        %.0f%%  (best published 11%%)"
          % (100.0 * pv.count("true") / len(pv)) if pv else "n/a")
    print("wrote FIG_sbdd_comparison.pdf / .png / .svg")

if __name__ == "__main__":
    main()
