from __future__ import annotations

import csv
import glob
import io as _io
import os
import sys
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S
import verify_divergence as VD
from structure_vs_literature import (PHARMACOPHORE, WARHEAD, groups, best_reference,
                                     campaign, _f)

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
OUT = os.path.join(HERE, "FIG_structure_panel")

FIGURE_PHARMACOPHORE = [(n, s) for (n, s) in PHARMACOPHORE if n.split(" (")[0] not in (
    "pyridine", "imidazole", "pyrazole", "phenol", "halogenated arene", "nitrile")]

PHARM_HL = (0.36, 0.72, 0.45)
WARHEAD_HL = (0.95, 0.55, 0.15)
MOL_PX = (560, 330)
C_SIG = "#2a7d4f"

def draw(mol, highlights: Dict[int, Tuple[float, float, float]], px=MOL_PX):
    from rdkit.Chem import AllChem
    m = Chem.Mol(mol)
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(px[0], px[1])
    o = d.drawOptions()
    o.bondLineWidth = 2
    o.highlightBondWidthMultiplier = 12
    o.fixedFontSize = 13
    o.clearBackground = True
    idx = sorted(highlights)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=idx,
                                       highlightAtomColors=highlights)
    d.FinishDrawing()
    return mpimg.imread(_io.BytesIO(d.GetDrawingText()))

def matched_atoms(mol, table) -> Tuple[set, List[str]]:
    found, names = set(), []
    for name, sma in table:
        q = Chem.MolFromSmarts(sma)
        if q is None:
            continue
        hits = mol.GetSubstructMatches(q)
        if hits:
            names.append(name)
            for h in hits:
                found |= set(int(i) for i in h)
    return found, names

def warhead_atoms(mol) -> Tuple[set, List[str]]:
    found, names = set(), []
    for name, sma in WARHEAD:
        q = Chem.MolFromSmarts(sma)
        if q is None:
            continue
        hits = mol.GetSubstructMatches(q)
        if hits:
            names.append(name)
            for h in hits:
                found |= set(int(i) for i in h)
    return found, names

def mcs_atoms(lead, ref) -> Tuple[set, int]:
    try:
        res = rdFMCS.FindMCS([lead, ref], timeout=20, ringMatchesRingOnly=True,
                             completeRingsOnly=True)
    except Exception:
        return set(), 0
    if not res or res.numAtoms == 0:
        return set(), 0
    q = Chem.MolFromSmarts(res.smartsString)
    if q is None:
        return set(), res.numAtoms
    m = lead.GetSubstructMatch(q)
    return set(int(i) for i in m), res.numAtoms

def active_pool() -> Dict[str, List[str]]:
    import re
    master = os.path.join(HERE, "master_results.csv")
    out: Dict[str, List[str]] = {}
    if not os.path.exists(master):
        return out
    for r in csv.DictReader(open(master, newline="", encoding="utf-8")):
        if r.get("kind") not in ("reference", "drug") or not r.get("smiles"):
            continue
        case = re.sub(r"_panel$", "", str(r.get("case") or "").split("__")[0])
        out.setdefault(case, []).append(r["smiles"])
    return out

def decoy_null(case: str, ref, pool: Dict[str, List[str]], n: int = 60) -> List[float]:
    import random
    src = [s for c, v in pool.items() if c != case for s in v]
    random.Random(20240517).shuffle(src)
    vals = []
    for s in src[:n]:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        _, k = mcs_atoms(m, ref)
        vals.append(100.0 * k / ref.GetNumHeavyAtoms())
    return vals

def collect(limit: int = 8) -> List[Dict]:
    refs = best_reference()
    pool = active_pool()
    out = []
    for path in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(path).replace("panel_plan_", "").replace(".csv", "")
        rows = campaign(path)
        if not rows or case not in refs:
            continue
        final = rows[-1]
        lead = Chem.MolFromSmiles(final.get("smiles") or "")
        ref = Chem.MolFromSmiles(refs[case]["smiles"] or "")
        if lead is None or ref is None:
            continue
        core, n_mcs = mcs_atoms(lead, ref)
        wh_idx, wh_names = warhead_atoms(lead)
        _, ref_wh = warhead_atoms(ref)
        cov_round = next((r for r in rows
                          if str(r.get("goal_feature") or "").lower().startswith("coval")), None)
        lead_ph = groups(lead, FIGURE_PHARMACOPHORE)
        ref_ph = groups(ref, FIGURE_PHARMACOPHORE)

        pct = (100.0 * n_mcs / ref.GetNumHeavyAtoms()) if n_mcs else 0.0
        null = decoy_null(case, ref, pool)
        p = ((sum(1 for x in null if x >= pct) + 1) / float(len(null) + 1)) if null else None

        scaf = Chem.MolFromSmiles(rows[0].get("scaffold_smiles") or "")
        contained = bool(scaf is not None and lead.HasSubstructMatch(scaf))

        out.append({
            "pct": pct, "null": null, "p": p, "contained": contained,
            "scaffold_ha": (scaf.GetNumHeavyAtoms() if scaf is not None else None),
            "lead_ha": lead.GetNumHeavyAtoms(),
            "case": case, "lead": lead, "ref": ref, "ref_name": refs[case]["name"],
            "core": core, "n_mcs": n_mcs, "wh_idx": wh_idx, "wh_names": wh_names,
            "ref_wh": ref_wh, "cov_residue": (cov_round or {}).get("goal_residue"),
            "kept": sorted(set(ref_ph) & set(lead_ph)),
            "added": sorted(set(lead_ph) - set(ref_ph)),
            "lead_aff": _f(final.get("affinity_target")),
            "ref_aff": refs[case]["affinity"],
            "rounds": len(rows),
        })
    out.sort(key=lambda d: (not d["wh_names"], -d["n_mcs"]))
    return out[:limit]

def main(limit: int = 8) -> int:
    S.use()
    items = collect(limit)
    if not items:
        print("no case has both a campaign and a published compound")
        return 1

    n = len(items)
    fig = plt.figure(figsize=(S.COL_2, 1.62 * n), constrained_layout=False)
    gs = fig.add_gridspec(n, 3, width_ratios=[1, 1.15, 0.72], hspace=0.46, wspace=0.05,
                          left=0.015, right=0.985, top=0.962, bottom=0.022)

    for i, d in enumerate(items):
        axl = fig.add_subplot(gs[i, 0])
        axr = fig.add_subplot(gs[i, 1])
        axe = fig.add_subplot(gs[i, 2])
        for ax in (axl, axr):
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)

        ref_ph_idx, _ = matched_atoms(d["ref"], FIGURE_PHARMACOPHORE)
        ref_wh_idx, _ = matched_atoms(d["ref"], WARHEAD)
        hl_ref = {a: PHARM_HL for a in ref_ph_idx}
        hl_ref.update({a: WARHEAD_HL for a in ref_wh_idx})
        axl.imshow(draw(d["ref"], hl_ref))
        axl.set_title(S.clean("%s   published: %s" % (d["case"], d["ref_name"][:26])),
                      fontsize=6.5, loc="left", pad=2, color=S.C_TEXT)

        lead_ph_idx, _ = matched_atoms(d["lead"], FIGURE_PHARMACOPHORE)
        hl = {a: PHARM_HL for a in lead_ph_idx}
        hl.update({a: WARHEAD_HL for a in d["wh_idx"]})
        axr.imshow(draw(d["lead"], hl))
        pct = d["pct"]
        axr.set_title(S.clean("generated lead, %d rounds   shared core %.0f%% of the published "
                              "compound" % (d["rounds"], pct)),
                      fontsize=6.5, loc="left", pad=2, color=S.C_TEXT)

        bits = []
        if d["kept"]:
            bits.append("kept " + ", ".join(x.split(" (")[0] for x in d["kept"][:3]))
        if d["wh_names"]:
            res = (" aimed at %s" % d["cov_residue"]) if d["cov_residue"] else ""
            if d["ref_wh"]:
                bits.append("covalent %s%s; the published compound uses %s at the same residue"
                            % (d["wh_names"][0].split(" (")[0], res,
                               d["ref_wh"][0].split(" (")[0]))
            else:
                bits.append("covalent %s%s, which the published compound does not have"
                            % (d["wh_names"][0].split(" (")[0], res))
        elif d["added"]:
            bits.append("added " + ", ".join(x.split(" (")[0] for x in d["added"][:2]))
        axr.text(0.0, -0.08, S.clean(".  ".join(bits)), transform=axr.transAxes,
                 fontsize=5.8, va="top", ha="left", color="#333333")
        axl.text(0.0, -0.08,
                 S.clean("docking %.2f vs %.2f kcal/mol, context only"
                         % (d["ref_aff"], d["lead_aff"]) if d["lead_aff"] is not None
                         else "docking %.2f kcal/mol" % d["ref_aff"]),
                 transform=axl.transAxes, fontsize=5.8, va="top", ha="left", color="#777777")

        axe.set_xlim(0, 100)
        axe.set_ylim(-2.35, 1.78)
        axe.set_xticks([]); axe.set_yticks([])
        for sp in axe.spines.values():
            sp.set_visible(False)

        axe.text(0, 1.72, S.clean("shared core against chance"), fontsize=5.4,
                 color=S.C_TEXT, fontweight="bold", va="top")
        axe.plot([0, 100], [0.42, 0.42], linewidth=0.6, color="#bbbbbb", zorder=1)
        for v in d["null"]:
            axe.plot([v, v], [0.26, 0.58], linewidth=0.45, color="#c9d2da", zorder=2)
        if d["null"]:
            med = sorted(d["null"])[len(d["null"]) // 2]
            axe.plot([med, med], [0.18, 0.66], linewidth=1.0, color="#8c8c8c", zorder=3)
            axe.text(med, -0.02, S.clean("null %.0f%%" % med), fontsize=4.3, ha="center",
                     va="top", color="#8c8c8c")
        axe.plot([pct, pct], [0.10, 0.76], linewidth=1.7,
                 color=(C_SIG if (d["p"] is not None and d["p"] < 0.05) else "#1b6ca8"), zorder=4)
        axe.text(pct, 0.86, S.clean("%.0f%%" % pct), fontsize=5.6, ha="center", va="bottom",
                 color=S.C_TEXT, fontweight="bold")
        if d["p"] is not None:
            axe.text(0, -0.52, S.clean("p = %.3f against %d decoys%s"
                                       % (d["p"], len(d["null"]),
                                          ", significant" if d["p"] < 0.05 else "")),
                     fontsize=4.6, va="top",
                     color=(C_SIG if d["p"] < 0.05 else "#777777"))

        marks = []
        if d["contained"]:
            marks.append("scaffold intact, %d to %d heavy atoms"
                         % (d["scaffold_ha"] or 0, d["lead_ha"]))
        else:
            marks.append("scaffold NOT intact")
        if d["wh_names"] and d["cov_residue"]:
            _num = "".join(c for c in str(d["cov_residue"]) if c.isdigit())
            _chk = VD.check(d["case"], int(_num)) if _num else {}
            if _chk.get("ok"):
                marks.append("warhead at %s, %s"
                             % (d["cov_residue"],
                                ("%s in the anti-target" % VD.three(_chk["anti_residue"]))
                                if _chk.get("divergent") else "CONSERVED in the anti-target"))
            else:
                marks.append("warhead at %s" % d["cov_residue"])
        axe.text(0, -1.02, S.clean("\n".join(marks)), fontsize=4.6, va="top",
                 color=(C_SIG if d["contained"] else "#b4472e"), linespacing=1.55)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), dpi=600 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf and %s.png  (%d cases)"
          % (os.path.relpath(OUT, HERE), os.path.relpath(OUT, HERE), n))
    for d in items:
        print("  %-16s shared %3d atoms   warhead %-22s rounds %d"
              % (d["case"][:16], d["n_mcs"],
                 (d["wh_names"][0].split(" (")[0] if d["wh_names"] else "none"), d["rounds"]))
    return 0

if __name__ == "__main__":
    raise SystemExit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
