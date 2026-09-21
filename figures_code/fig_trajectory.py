from __future__ import annotations

import csv
import io as _io
import json
import os
import sys
import textwrap
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
OUT = os.path.join(HERE, "FIG_trajectory")

CASES = ["ABL1_SRC", "CA2_CA1", "EGFR_T790M", "KRAS_G12C"]

TITLES = {
    "ABL1_SRC": "ABL1 over SRC",
    "CA2_CA1": "Carbonic anhydrase II over carbonic anhydrase I",
    "EGFR_T790M": "EGFR T790M over wild-type EGFR",
    "KRAS_G12C": "KRAS G12C over wild-type KRAS",
}

NEW_HL = (0.95, 0.55, 0.15)
C_OK = "#2a7d4f"
C_GREY = "#8c8c8c"
MOL_PX = (470, 260)

def reasons_for(case: str) -> Dict[str, str]:
    path = os.path.join(RERUN, "run_%s.trace.json" % case)
    if not os.path.exists(path):
        return {}
    try:
        d = json.load(open(path, encoding="utf-8"))
    except Exception:
        return {}
    out: Dict[str, str] = {}
    for e in d.get("events") or []:
        if e.get("kind") != "reasoning":
            continue
        title = str(e.get("title") or "")
        if not title.startswith("aim at"):
            continue
        parts = title.split()
        if len(parts) >= 3:
            out[parts[2]] = str(e.get("detail") or "")
    return out

def rounds(case: str) -> List[Dict]:
    path = os.path.join(RERUN, "panel_plan_%s.csv" % case)
    if not os.path.exists(path):
        return []
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8"))
            if (r.get("smiles") or "").strip()]
    rows.sort(key=lambda r: int(r.get("round") or 0))
    return rows

def new_atoms(mol, parent) -> set:
    if parent is None or mol is None:
        return set()
    m = mol.GetSubstructMatch(parent)
    if not m:
        return set()
    return set(range(mol.GetNumAtoms())) - set(int(i) for i in m)

def draw(mol, hl: set, px=MOL_PX):
    m = Chem.Mol(mol)
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(px[0], px[1])
    o = d.drawOptions()
    o.bondLineWidth = 2
    o.fixedFontSize = 12
    o.clearBackground = True
    idx = sorted(hl)
    rdMolDraw2D.PrepareAndDrawMolecule(
        d, m, highlightAtoms=idx, highlightAtomColors={a: NEW_HL for a in idx})
    d.FinishDrawing()
    return _crop(mpimg.imread(_io.BytesIO(d.GetDrawingText())))

def _crop(img, pad: int = 6):
    try:
        a = img
        if a.ndim != 3:
            return img
        if a.shape[2] == 4:
            mask = a[:, :, 3] > 0.02
        else:
            mask = a[:, :, :3].min(axis=2) < 0.97
        ys, xs = np.where(mask)
        if ys.size == 0 or xs.size == 0:
            return img
        y0 = max(int(ys.min()) - pad, 0)
        y1 = min(int(ys.max()) + pad + 1, a.shape[0])
        x0 = max(int(xs.min()) - pad, 0)
        x1 = min(int(xs.max()) + pad + 1, a.shape[1])
        return a[y0:y1, x0:x1]
    except Exception:
        return img

def build_case(case: str, rws, reasons) -> str:
    scaf_smi = (rws[0].get("scaffold_smiles") or "").strip()
    seq = []
    if scaf_smi and Chem.MolFromSmiles(scaf_smi) is not None:
        seq.append({"smiles": scaf_smi, "round": "0", "_scaffold_row": True})
    seq.extend(rws)
    n_rounds = len(rws)
    rws = seq
    n = len(rws)
    fig = plt.figure(figsize=(S.COL_2, 0.62 * n + 0.62))
    gs = fig.add_gridspec(n + 1, 2, width_ratios=[1.0, 1.42], hspace=0.18, wspace=0.03,
                          left=0.012, right=0.988, top=0.975, bottom=0.015)

    axh = fig.add_subplot(gs[0, :])
    axh.axis("off"); axh.set_xlim(0, 1); axh.set_ylim(0, 1)
    axh.text(0.0, 0.55, S.clean(TITLES.get(case, case)), fontsize=8.6,
             fontweight="bold", color=S.C_TEXT, va="center")
    first = Chem.MolFromSmiles(rws[0].get("scaffold_smiles") or rws[0].get("smiles") or "")
    last = Chem.MolFromSmiles(rws[-1].get("smiles") or "")
    if first is not None and last is not None:
        axh.text(0.32, 0.55,
                 S.clean("%d rounds from a %d heavy-atom scaffold to %d, scaffold %s in the "
                         "final lead"
                         % (n_rounds, first.GetNumHeavyAtoms(), last.GetNumHeavyAtoms(),
                            "intact" if last.HasSubstructMatch(first) else "LOST")),
                 fontsize=6.0, color=C_GREY, va="center")
    axh.plot([0, 1], [0.06, 0.06], linewidth=0.7, color="#444444")

    parent = None
    for i, r in enumerate(rws):
        mol = Chem.MolFromSmiles(r.get("smiles") or "")
        if mol is None:
            continue
        axm = fig.add_subplot(gs[i + 1, 0])
        axt = fig.add_subplot(gs[i + 1, 1])
        for ax in (axm, axt):
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)

        hl = new_atoms(mol, parent)
        axm.imshow(draw(mol, hl))

        axt.set_xlim(0, 1); axt.set_ylim(0, 1); axt.axis("off")
        rnum = r.get("round") or "0"
        res = (r.get("goal_residue") or "").strip()
        feat = (r.get("goal_feature") or "").strip()
        head_t = ("round 0: scaffold, mined from the target's known actives"
                  if r.get("_scaffold_row")
                  else ("round %s: add %s aimed at %s" % (rnum, feat or "a group", res)
                        if res else "round %s" % rnum))
        axt.text(0.0, 0.94, S.clean(head_t), fontsize=6.4, fontweight="bold",
                 color=S.C_TEXT, va="top")

        contained = bool(parent is None or mol.HasSubstructMatch(parent))
        if r.get("_scaffold_row"):
            body = "%d heavy atoms. This is what every later round must keep."                 % mol.GetNumHeavyAtoms()
        else:
            body = ("%d heavy atoms, %d installed this round.  %s"
                    % (mol.GetNumHeavyAtoms(), len(hl),
                       "contains the previous round" if contained
                       else "DOES NOT CONTAIN THE PREVIOUS ROUND"))
        axt.text(0.0, 0.72, S.clean(body), fontsize=5.5,
                 color=(C_OK if contained else "#b4472e"), va="top")

        axis = (r.get("goal_axis") or "").strip()
        gate = " ".join((r.get("sel_gate_reason") or "").split())
        bits = []
        if axis:
            nice = {"agent_divergent": "position DIFFERS in the anti-target",
                    "agent_conserved": "position is CONSERVED in the anti-target"}.get(axis, axis)
            bits.append("engine logged this as: %s" % nice)
        if gate:
            bits.append("counter-screen, reported not blocking: %s" % gate[:170])
        if bits:
            axt.text(0.0, 0.53, S.clean(textwrap.fill("  |  ".join(bits), 94)), fontsize=5.0,
                     color=(C_OK if axis == "agent_divergent" else C_GREY), va="top",
                     linespacing=1.45)
        why = reasons.get(res) if res else None
        if why:
            axt.text(0.0, 0.18,
                     S.clean(textwrap.fill("planner: " + " ".join(why.split())[:300], 96)),
                     fontsize=4.9, color="#333333", va="top", linespacing=1.45)
        parent = mol

    path = "%s_%s" % (OUT, case)
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (path, ext), dpi=500 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return path

def main() -> int:
    S.use()
    made = 0
    for case in CASES:
        rws = rounds(case)
        if not rws:
            print("  %s: no campaign record" % case)
            continue
        path = build_case(case, rws, reasons_for(case))
        made += 1
        print("wrote %s.pdf and .png  (%d rounds)" % (os.path.relpath(path, HERE), len(rws)))
    if not made:
        print("no campaign records found")
        return 1
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
