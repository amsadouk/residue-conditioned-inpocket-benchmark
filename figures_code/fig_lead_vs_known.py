from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
MASTER = os.path.join(os.path.dirname(HERE), "data", "tables", "master_results.csv")
OUT = os.path.join(HERE, "figures_leads")

PANEL_W, PANEL_H = 520, 380
ADDED_RGB = (0.20, 0.65, 0.35)
CORE_RGB = (0.55, 0.60, 0.68)

def _f(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() in ("", "None", "nan", "-"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

def best_known(case: str) -> Optional[Dict]:
    if not os.path.exists(MASTER):
        return None
    best = None
    base = case.split("__")[0]
    for r in csv.DictReader(open(MASTER, newline="", encoding="utf-8")):
        if str(r.get("case") or "").split("__")[0] != base:
            continue
        if r.get("kind") not in ("reference", "drug"):
            continue
        a = _f(r.get("affinity_target"))
        if a is None:
            continue
        if best is None or a < best["affinity"]:
            best = {"affinity": a, "smiles": r.get("smiles"),
                    "name": (r.get("resolved_name") or r.get("name") or "reference"),
                    "le": _f(r.get("le_target"))}
    return best

def campaign(case_csv: str) -> List[Dict]:
    rows = []
    for r in csv.DictReader(open(case_csv, newline="", encoding="utf-8")):
        rows.append({
            "round": r.get("round"), "smiles": r.get("smiles"),
            "scaffold": r.get("scaffold_smiles"),
            "affinity": _f(r.get("affinity_target")),
            "margin": _f(r.get("margin")) if _f(r.get("margin")) is not None
                      else _f(r.get("counter_margin")),
            "residue": r.get("goal_residue"), "feature": r.get("goal_feature"),
            "divergent": _f(r.get("divergent_contacts")),
            "delta": _f(r.get("delta_this_round")),
            "heavy": r.get("heavy_atoms"),
        })
    rows = [r for r in rows if r["smiles"]]
    rows.sort(key=lambda r: int(r["round"] or 0))
    return rows

def added_atoms(lead_smiles: str, scaffold_smiles: str):
    lead = Chem.MolFromSmiles(lead_smiles or "")
    scaf = Chem.MolFromSmiles(scaffold_smiles or "")
    if lead is None:
        return None, set(), set()
    if scaf is None:
        return lead, set(), set()
    match = lead.GetSubstructMatch(scaf)
    if not match:
        return lead, set(), set()
    core = set(int(i) for i in match)
    added = set(a.GetIdx() for a in lead.GetAtoms()) - core
    return lead, core, added

def draw_pair(case: str, known: Optional[Dict], rounds: List[Dict], path: str) -> Optional[str]:
    if not rounds:
        return None
    lead_row = min((r for r in rounds if r["affinity"] is not None),
                   key=lambda r: r["affinity"], default=None)
    if lead_row is None:
        return None
    lead, core, added = added_atoms(lead_row["smiles"], lead_row["scaffold"])
    if lead is None:
        return None

    mols, legends, highlights, hcolors = [], [], [], []
    if known and known.get("smiles"):
        km = Chem.MolFromSmiles(known["smiles"])
        if km is not None:
            mols.append(km)
            legends.append("%s\nbest known: %.2f kcal/mol"
                           % (str(known["name"])[:34], known["affinity"]))
            highlights.append([])
            hcolors.append({})
    mols.append(lead)
    legends.append("generated lead (round %s)\n%.2f kcal/mol   margin %s"
                   % (lead_row["round"], lead_row["affinity"],
                      ("%+.2f" % lead_row["margin"]) if lead_row["margin"] is not None else "n/a"))
    highlights.append(sorted(core | added))
    hcolors.append({**{i: CORE_RGB for i in core}, **{i: ADDED_RGB for i in added}})

    for m in mols:
        AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(mols, molsPerRow=len(mols),
                               subImgSize=(PANEL_W, PANEL_H), legends=legends,
                               highlightAtomLists=highlights,
                               highlightAtomColors=hcolors, useSVG=True)
    svg = img.data if hasattr(img, "data") else str(img)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(svg)
    return path

def rationale(case: str, rounds: List[Dict], known: Optional[Dict]) -> str:
    out = ["%s" % case]
    if known:
        out.append("  best known compound: %s at %.2f kcal/mol" % (known["name"], known["affinity"]))
    out.append("  scaffold: %s" % (rounds[0]["scaffold"] if rounds else "?"))
    for r in rounds:
        cov = " [COVALENT: no scoring function here forms a bond, so this margin does not " \
              "measure the selectivity the warhead was chosen for]" \
              if str(r.get("feature") or "").lower().startswith("coval") else ""
        out.append("  r%-2s aimed %-8s %-14s  affinity %-8s margin %-7s  %s%s"
                   % (r["round"], str(r["residue"] or "?"), str(r["feature"] or "?"),
                      ("%.2f" % r["affinity"]) if r["affinity"] is not None else "n/a",
                      ("%+.2f" % r["margin"]) if r["margin"] is not None else "n/a",
                      ("%+.2f this round" % r["delta"]) if r["delta"] is not None else "",
                      cov))
    return "\n".join(out)

def main(limit: int = 8) -> int:
    import glob
    os.makedirs(OUT, exist_ok=True)
    cases = sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv")))
    made = 0
    notes = []
    for path in cases:
        if made >= limit:
            break
        case = os.path.basename(path).replace("panel_plan_", "").replace(".csv", "")
        rounds = campaign(path)
        if not rounds:
            continue
        known = best_known(case)
        svg = os.path.join(OUT, "lead_vs_known_%s.svg" % case)
        if draw_pair(case, known, rounds, svg):
            made += 1
            notes.append(rationale(case, rounds, known))
            print("wrote %s" % os.path.relpath(svg, HERE))
    if notes:
        with open(os.path.join(OUT, "RATIONALE.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n\n".join(notes) + "\n")
        print("\n" + "\n\n".join(notes))
    print("\n%d panel(s) in %s" % (made, os.path.relpath(OUT, HERE)))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
