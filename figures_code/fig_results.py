import os, csv, json, io, traceback
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from scipy import stats as sps

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem.Draw import rdMolDraw2D

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, GRID = "#1f2328", "#6a737d", "#DDE1E5"
GEN, REF = "#B5651D", "#3D5A80"
GOOD, BAD = "#2F6F62", "#8B2E2E"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                     "axes.linewidth": 0.7, "axes.edgecolor": "#B9BEC4"})

AXES = [("selectivity_ratio", "pocket specificity"),
        ("margin", "selectivity margin (kcal/mol)"),
        ("affinity_target", "docking score (kcal/mol)"),
        ("le_target", "ligand efficiency")]
LOWER_IS_BETTER = {"affinity_target"}

def f(x):
    try:
        if x is None or str(x).strip() in ("", "None", "nan"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

def load():
    rows = list(csv.DictReader(open(os.path.join(HERE, "FINAL_benchmark.csv"),
                                    newline="", encoding="utf-8")))

    def live(r):
        return (str(r.get("superseded_by_redock") or "").lower() != "true"
                and str(r.get("is_peptide_like") or "").lower() != "true")
    return [r for r in rows if live(r)]

def short(case):
    return case.replace("__selgoal", "").replace("__anchor", " (anchor)")

def draw_mol(smi, w=330, h=250, highlight=None):
    m = Chem.MolFromSmiles(smi or "")
    if m is None:
        return None
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    o = d.drawOptions()
    o.clearBackground = True
    o.bondLineWidth = 1.7
    o.fixedFontSize = 15
    o.minFontSize = 11
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")

def fig2(rows):
    S = json.load(open(os.path.join(HERE, "stats_multi.json"), encoding="utf-8"))
    cases = [c for c in S if "__anchor" not in c]
    cases.sort(key=lambda c: -(S[c].get("selectivity_ratio", {}).get("rb") or -9))
    if "PTGS2_PTGS1__anchor" in S:
        cases.append("PTGS2_PTGS1__anchor")

    fig, axs = plt.subplots(1, 4, figsize=(15.2, 0.46 * len(cases) + 2.5), dpi=300,
                            sharey=True)
    fig.subplots_adjust(left=0.155, right=0.985, top=0.80, bottom=0.12, wspace=0.13)
    ypos = np.arange(len(cases))[::-1]

    for j, (key, title) in enumerate(AXES):
        ax = axs[j]
        ax.axvline(0, color="#9AA0A6", lw=0.9, zorder=1)
        ax.axvspan(0, 1.06, color=GEN, alpha=0.045, zorder=0)
        ax.axvspan(-1.06, 0, color=REF, alpha=0.045, zorder=0)
        for y, c in zip(ypos, cases):
            st = S.get(c, {}).get(key) or {}
            rb, p = st.get("rb"), st.get("p")
            if rb is None:
                ax.text(0, y, "not testable", fontsize=6.6, color=MUTED, ha="center",
                        va="center", style="italic")
                continue
            val = -rb if key in LOWER_IS_BETTER else rb
            sig = p is not None and p < 0.05
            col = GEN if val > 0 else REF
            ax.plot([0, val], [y, y], color=col, lw=2.6 if sig else 1.4,
                    alpha=1.0 if sig else 0.5, solid_capstyle="round", zorder=3)
            ax.plot([val], [y], marker="o", ms=7.5 if sig else 5.2, color=col,
                    markeredgecolor="white", markeredgewidth=1.0,
                    alpha=1.0 if sig else 0.55, zorder=4)
            if sig:
                ax.text(val + (0.055 if val > 0 else -0.055), y,
                        "p=%.4f" % p if p >= 0.0001 else "p<0.0001",
                        fontsize=6.9, fontweight="bold", color=col,
                        ha="left" if val > 0 else "right", va="center", zorder=5)
        ax.set_xlim(-1.32, 1.32)
        ax.set_ylim(-0.75, len(cases) - 0.25)
        ax.set_yticks(ypos)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_xticklabels(["-1", "", "0", "", "+1"], fontsize=7.4)
        ax.set_title(title, fontsize=9.4, fontweight="bold", color=INK, pad=9)
        ax.tick_params(length=2.5, colors=MUTED)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)

    axs[0].set_yticklabels([short(c) for c in cases], fontsize=8.4, color=INK)
    for t, c in zip(axs[0].get_yticklabels(), cases):
        if "anchor" in c:
            t.set_color(MUTED)
            t.set_style("italic")

    fig.text(0.155, 0.945, "Generated molecules versus each target's own known compounds",
             fontsize=15, fontweight="bold", color=INK)
    fig.text(0.155, 0.912,
             "Rank-biserial correlation, Mann-Whitney U, one-sided. Right of zero favours the "
             "generated molecules; +1.00 is complete separation.\n"
             "Bold markers are p < 0.05. n = 3 generated against 10 references per case. "
             "The italic row is the conserved-anchor control arm.",
             fontsize=8.6, color=MUTED, va="top")
    fig.text(0.155, 0.045, "generated better  →", fontsize=8.2, color=GEN,
             fontweight="bold")
    fig.text(0.985, 0.045, "←  reference better", fontsize=8.2, color=REF,
             fontweight="bold", ha="right")
    p = os.path.join(HERE, "Figure2_benchmark.png")
    fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    return p

def fig3(rows):
    leads = [r for r in rows if r.get("kind") == "generated"
             and "__selgoal" in (r.get("case") or "") and (r.get("smiles") or "").strip()]
    by = {}
    for r in leads:
        by.setdefault(r["case"], []).append(r)
    for c in by:
        by[c].sort(key=lambda r: int(f(r.get("round")) or 0))
    cases = sorted(by, key=lambda c: -max((f(x.get("selectivity_ratio")) or 0)
                                          for x in by[c]))
    ncol = max(len(v) for v in by.values())

    fig = plt.figure(figsize=(3.55 * ncol + 1.9, 2.52 * len(cases) + 1.7), dpi=300)
    fig.subplots_adjust(left=0.085, right=0.99, top=0.915, bottom=0.022,
                        hspace=0.30, wspace=0.10)

    for i, c in enumerate(cases):
        for j in range(ncol):
            ax = fig.add_subplot(len(cases), ncol, i * ncol + j + 1)
            ax.axis("off")
            if j >= len(by[c]):
                continue
            r = by[c][j]
            img = draw_mol(r["smiles"], 360, 260)
            if img is not None:
                ax.imshow(img, interpolation="lanczos")
            div = (r.get("goal_axis") or "").startswith("divergent")
            col = GOOD if div else MUTED
            ax.text(0.5, 1.10, "round %s   ·   %s / %s"
                    % (r.get("round") or "?", r.get("goal_residue") or "?",
                       r.get("goal_feature") or "?"),
                    transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=col,
                    ha="center", va="bottom")
            ax.text(0.5, 1.035,
                    "divergent target" if div else "conserved fallback",
                    transform=ax.transAxes, fontsize=7.0, color=col, ha="center",
                    va="bottom", style="italic")
            spec = r.get("selectivity_ratio")
            ax.text(0.5, -0.045,
                    "%s kcal/mol    margin %s    spec %s    LE %s"
                    % (r.get("affinity_target") or "-", r.get("margin") or "-",
                       spec if spec not in (None, "", "None") else "-",
                       r.get("le_target") or "-"),
                    transform=ax.transAxes, fontsize=7.6, color=INK, ha="center", va="top")
            if j == 0:
                ax.text(-0.10, 0.5, short(c), transform=ax.transAxes, fontsize=10.6,
                        fontweight="bold", color=INK, rotation=90, ha="center", va="center")

    fig.text(0.085, 0.975, "Every molecule the loop carried forward",
             fontsize=16, fontweight="bold", color=INK)
    fig.text(0.085, 0.951,
             "Rows are targets, columns are growth rounds. Each structure is the survivor the "
             "engine selected from that round's Pareto front and docked into the pinned pocket.",
             fontsize=9.0, color=MUTED, va="top")
    p = os.path.join(HERE, "Figure3_molecules.png")
    fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    return p

def fig4(rows):
    gen = [r for r in rows if r.get("case") == "PTGS2_PTGS1__selgoal"
           and r.get("kind") == "generated"]
    ref = [r for r in rows if r.get("case") == "PTGS2_PTGS1__selgoal"
           and r.get("kind") == "reference" and f(r.get("affinity_target")) is not None]
    seen, uniq = set(), []
    for r in sorted(gen, key=lambda r: int(f(r.get("round")) or 0)):
        k = Chem.CanonSmiles(r["smiles"]) if r.get("smiles") else None
        if k and k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    from rdkit.Chem import AllChem, DataStructs
    lead = Chem.MolFromSmiles(uniq[0]["smiles"]) if uniq else None
    if lead is not None:
        fpl = AllChem.GetMorganFingerprintAsBitVect(lead, 2, nBits=2048)

        def sim(r):
            m = Chem.MolFromSmiles(r.get("smiles") or "")
            if m is None:
                return -1.0
            return DataStructs.TanimotoSimilarity(
                fpl, AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048))
        for r in ref:
            r["_sim"] = sim(r)
        ref.sort(key=lambda r: -r["_sim"])
    else:
        ref.sort(key=lambda r: -(f(r.get("pchembl")) or 0))
    show = uniq + ref[:3]

    fig = plt.figure(figsize=(4.0 * len(show) + 1.0, 6.0), dpi=300)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.735, bottom=0.20, wspace=0.08)
    for j, r in enumerate(show):
        ax = fig.add_subplot(1, len(show), j + 1)
        ax.axis("off")
        img = draw_mol(r["smiles"], 430, 320)
        if img is not None:
            ax.imshow(img, interpolation="lanczos")
        isgen = r.get("kind") == "generated"
        col = GEN if isgen else REF
        head = ("GENERATED  round %s" % r.get("round")) if isgen else "KNOWN COMPOUND"
        ax.text(0.5, 1.10, head, transform=ax.transAxes, fontsize=9.6,
                fontweight="bold", color=col, ha="center", va="bottom")
        if isgen:
            sub = "aimed %s / %s" % (r.get("goal_residue"), r.get("goal_feature"))
        else:
            nm = (r.get("resolved_name") or r.get("label") or "").strip()
            if len(nm) > 30:
                cut = nm[:30].rsplit(" ", 1)[0] if " " in nm[:30] else nm[:30]
                nm = cut + "…"
            s = r.get("_sim")
            sub = nm + (("\nTanimoto %.2f to the lead" % s) if isinstance(s, float)
                        and s >= 0 else "")
        ax.text(0.5, 1.045, sub, transform=ax.transAxes, fontsize=7.6, color=col,
                ha="center", va="bottom", style="italic")
        pc = r.get("pchembl")
        ax.text(0.5, -0.04,
                "%s kcal/mol\nmargin %s   spec %s   LE %s%s"
                % (r.get("affinity_target") or "-", r.get("margin") or "-",
                   r.get("selectivity_ratio") or "-", r.get("le_target") or "-",
                   ("\npChEMBL %s" % pc) if pc not in (None, "", "None") else ""),
                transform=ax.transAxes, fontsize=8.0, color=INK, ha="center", va="top")

    fig.text(0.02, 0.965, "COX-2 over COX-1: the engine finds the chemotype and misses the "
             "pharmacophore", fontsize=15.5, fontweight="bold", color=INK)
    fig.text(0.02, 0.905,
             "Both arms share the 1,2-diarylpyrrole core of the coxibs. The known compound "
             "carries the para-sulfonamide that reaches Arg513 in the COX-2 side pocket and "
             "drives isoform selectivity;\nthe generated molecules replace it with plain aryl "
             "and add a cyclopropyl elsewhere. That single difference explains both results: "
             "higher ligand efficiency because it is smaller,\nlower selectivity margin "
             "because the group that buys COX-2 over COX-1 is absent. Rounds 1 and 3 are the "
             "same molecule - round 3 failed to improve and the front collapsed back.",
             fontsize=8.5, color=MUTED, va="top")
    p = os.path.join(HERE, "Figure4_cox2.png")
    fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    return p

def fig5(rows):
    pts = {}
    for r in rows:
        if r.get("kind") not in ("reference", "drug"):
            continue
        s, pc = f(r.get("affinity_target")), f(r.get("pchembl"))
        if s is None or pc is None:
            continue
        gene = (r.get("gene") or "").strip() or \
            (r.get("case") or "").replace("__selgoal", "").replace("__anchor", "").split("_")[0]
        key = "%s %s" % (gene or "?", (r.get("target_pdb") or "?"))
        ident = (r.get("inchikey") or "").strip() or Chem.CanonSmiles(r["smiles"]) \
            if r.get("smiles") else (r.get("label") or "")
        pts.setdefault(key, {})[ident] = (s, pc)
    keys = sorted(pts, key=lambda k: -len(pts[k]))
    keys = [k for k in keys if len(pts[k]) >= 4]

    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 3.15 * nrow + 1.35), dpi=300)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.845, bottom=0.085,
                        hspace=0.46, wspace=0.30)
    axs = np.atleast_1d(axs).ravel()
    nbad = 0
    for ax, k in zip(axs, keys):
        v = list(pts[k].values())
        x = np.array([a for a, _ in v])
        y = np.array([b for _, b in v])
        rho, p = sps.spearmanr(x, y)
        backwards = rho > 0
        nbad += 1 if backwards else 0
        col = BAD if backwards else GOOD
        ax.scatter(x, y, s=46, facecolor=col, edgecolor="white", linewidth=0.9,
                   alpha=0.92, zorder=3)
        if len(x) > 2:
            b, a = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 10)
            ax.plot(xs, a + b * xs, color=col, lw=1.5, ls="--", alpha=0.75, zorder=2)
        ax.set_title(k, fontsize=9.6, fontweight="bold", color=INK, pad=6)
        ax.text(0.5, 1.005, "Spearman rho = %+.3f    p = %.3f    n = %d" % (rho, p, len(x)),
                transform=ax.transAxes, fontsize=7.6, color=col, ha="center", va="bottom",
                fontweight="bold" if backwards else "normal")
        if backwards:
            ax.text(0.97, 0.05, "BACKWARDS", transform=ax.transAxes, fontsize=8.6,
                    color=BAD, fontweight="bold", ha="right", va="bottom")
        ax.set_xlabel("docking score (kcal/mol)", fontsize=8.0, color=MUTED)
        ax.set_ylabel("measured pChEMBL", fontsize=8.0, color=MUTED)
        ax.tick_params(labelsize=7.4, length=2.5, colors=MUTED)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(color=GRID, lw=0.5, zorder=0)
    for ax in axs[len(keys):]:
        ax.axis("off")

    fig.text(0.075, 0.965, "Does the score track real potency? Per receptor, never pooled.",
             fontsize=15.5, fontweight="bold", color=INK)
    fig.text(0.075, 0.905,
             "A more negative docking score should go with a HIGHER measured pChEMBL, so a "
             "NEGATIVE slope means the oracle is working. %d of %d receptors run the wrong "
             "way.\nPooling these would average opposite-signed relationships measured on "
             "different score scales and report a mild nothing; it is not done here. "
             "Where the slope is inverted,\nno ranking on that target is supported by this "
             "benchmark however the per-case comparison comes out."
             % (nbad, len(keys)),
             fontsize=8.5, color=MUTED, va="top")
    p = os.path.join(HERE, "Figure5_oracle.png")
    fig.savefig(p, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.22)
    plt.close(fig)
    return p

if __name__ == "__main__":
    rows = load()
    print("live rows: %d" % len(rows))
    for name, fn in (("Fig2", fig2), ("Fig3", fig3), ("Fig4", fig4), ("Fig5", fig5)):
        try:
            print("%s -> %s" % (name, os.path.basename(fn(rows))))
        except Exception:
            print("%s FAILED" % name)
            traceback.print_exc()
