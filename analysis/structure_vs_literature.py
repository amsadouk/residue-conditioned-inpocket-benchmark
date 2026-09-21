from __future__ import annotations

import csv
import glob
import os
import sys
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import Descriptors, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
MASTER = os.path.join(os.path.dirname(HERE), "data", "tables", "master_results.csv")
OUT = os.path.join(HERE, "STRUCTURE_vs_literature.txt")

ROLES = {
    "primary sulfonamide (Zn binder)": {
        "roles": ["primary_recognition_anchor", "metal_coordination"],
        "primary_recognition_anchor": True},
    "aminopyrimidine (hinge)": {
        "roles": ["primary_recognition_anchor", "hinge_binding"],
        "primary_recognition_anchor": True},
    "quinazoline (hinge)": {
        "roles": ["primary_recognition_anchor", "hinge_binding"],
        "primary_recognition_anchor": True},
    "aminoquinazoline (hinge)": {
        "roles": ["primary_recognition_anchor", "hinge_binding"],
        "primary_recognition_anchor": True},
    "anilide / secondary amide": {
        "roles": ["primary_recognition_anchor", "hydrogen_bonding"],
        "primary_recognition_anchor": True},
    "piperazine (solubilising)": {
        "roles": ["physicochemical_modulator", "peripheral_binding", "metabolic_soft_spot"],
        "primary_recognition_anchor": False, "literature_supported": True},
    "piperidine": {
        "roles": ["physicochemical_modulator", "peripheral_binding"],
        "primary_recognition_anchor": False},
    "morpholine": {
        "roles": ["physicochemical_modulator", "peripheral_binding"],
        "primary_recognition_anchor": False},
    "halogenated arene": {
        "roles": ["physicochemical_modulator", "peripheral_binding"],
        "primary_recognition_anchor": False},
}

def is_anchor(name: str) -> bool:
    return bool(ROLES.get(name, {}).get("primary_recognition_anchor", True))

def role_note(name: str) -> str:
    r = ROLES.get(name)
    if not r or r.get("primary_recognition_anchor", True):
        return ""
    return ", ".join(r.get("roles", []))

AUXILIARY = {k for k, v in ROLES.items() if not v.get("primary_recognition_anchor", True)}

PHARMACOPHORE: List[Tuple[str, str]] = [
    ("primary sulfonamide (Zn binder)", "[SX4](=O)(=O)[NX3;H2]"),
    ("secondary sulfonamide", "[SX4](=O)(=O)[NX3;H1]"),
    ("methylsulfone", "[SX4](=O)(=O)[CH3]"),
    ("hydroxamic acid", "C(=O)[NX3][OX2H1]"),
    ("carboxylic acid", "C(=O)[OX2H1]"),
    ("tetrazole (acid bioisostere)", "c1nn[nH]n1"),
    ("primary amide", "C(=O)[NX3;H2]"),
    ("anilide / secondary amide", "C(=O)[NX3;H1][c,C]"),
    ("urea", "[NX3][CX3](=O)[NX3]"),
    ("thiourea", "[NX3][CX3](=S)[NX3]"),
    ("aminopyrimidine (hinge)", "[NX3;H1]c1ncccn1"),
    ("quinazoline (hinge)", "c1ncnc2ccccc12"),
    ("aminoquinazoline (hinge)", "[NX3;H1]c1ncnc2ccccc12"),
    ("imidazole", "c1c[nH]cn1"),
    ("pyrazole", "c1cc[nH]n1"),
    ("pyridine", "c1ccncc1"),
    ("piperazine (solubilising)", "N1CCNCC1"),
    ("piperidine", "N1CCCCC1"),
    ("morpholine", "C1COCCN1"),
    ("oxime", "[CX3]=[NX2][OX2H1]"),
    ("phenol", "c[OX2H1]"),
    ("nitrile", "C#N"),
    ("halogenated arene", "c[F,Cl,Br,I]"),
]

WARHEAD: List[Tuple[str, str]] = [
    ("acrylamide (Michael acceptor)", "C=CC(=O)[NX3]"),
    ("propiolamide (alkyne)", "C#CC(=O)[NX3]"),
    ("vinyl sulfone", "C=C[SX4](=O)(=O)"),
    ("chloroacetamide", "[Cl]C[CX3](=O)[NX3]"),
    ("alpha-halo carbonyl", "[F,Cl,Br]C[CX3]=O"),
    ("epoxide", "C1OC1"),
    ("acrylonitrile (Michael acceptor)", "[CX3]=[CX3]C#N"),
]

def _f(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() in ("", "None", "nan", "-"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

def groups(mol, table) -> List[str]:
    out = []
    for name, sma in table:
        q = Chem.MolFromSmarts(sma)
        if q is not None and mol.HasSubstructMatch(q):
            out.append(name)
    return out

def shared_core(a, b) -> Dict:
    try:
        res = rdFMCS.FindMCS([a, b], timeout=20, ringMatchesRingOnly=True,
                             completeRingsOnly=True,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    except Exception:
        return {"atoms": 0, "smarts": None}
    if not res or res.canceled or res.numAtoms == 0:
        return {"atoms": 0, "smarts": None}
    return {"atoms": res.numAtoms, "bonds": res.numBonds, "smarts": res.smartsString}

def best_reference() -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    if not os.path.exists(MASTER):
        return out
    for r in csv.DictReader(open(MASTER, newline="", encoding="utf-8")):
        if r.get("kind") not in ("reference", "drug"):
            continue
        case = str(r.get("case") or "").split("__")[0]
        a = _f(r.get("affinity_target"))
        if a is None or not r.get("smiles"):
            continue
        if case not in out or a < out[case]["affinity"]:
            out[case] = {"affinity": a, "smiles": r.get("smiles"),
                         "name": (r.get("resolved_name") or r.get("name") or "reference")}
    return out

def campaign(path: str) -> List[Dict]:
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8")) if r.get("smiles")]
    rows.sort(key=lambda r: int(r.get("round") or 0))
    return rows

def report(case: str, rows: List[Dict], ref: Optional[Dict]) -> str:
    final = rows[-1]
    lead = Chem.MolFromSmiles(final["smiles"])
    if lead is None:
        return "%s: the final lead does not parse\n" % case
    out = []
    out.append("=" * 88)
    out.append(case)
    out.append("=" * 88)

    out.append("INSTALLATION LOG  (what each round aimed at)")
    for r in rows:
        out.append("  r%-2s  %-9s %-16s  affinity %-8s  margin %-7s  %s"
                   % (r.get("round"), r.get("goal_residue") or "?",
                      r.get("goal_feature") or "?",
                      r.get("affinity_target") or "n/a", r.get("margin") or "n/a",
                      "COVALENT" if str(r.get("goal_feature") or "").lower().startswith("coval")
                      else ""))
    out.append("")

    scaf = Chem.MolFromSmiles(rows[0].get("scaffold_smiles") or "")
    out.append("FINAL LEAD")
    out.append("  %s" % final["smiles"])
    out.append("  %d heavy atoms, MW %.0f, %d rings"
               % (lead.GetNumHeavyAtoms(), Descriptors.MolWt(lead),
                  lead.GetRingInfo().NumRings()))
    if scaf is not None:
        kept = lead.HasSubstructMatch(scaf)
        out.append("  starting scaffold %s (%d heavy atoms) is %s in the final lead"
                   % (Chem.MolToSmiles(scaf), scaf.GetNumHeavyAtoms(),
                      "STILL PRESENT WHOLE" if kept else "NOT PRESENT - the scheme was violated"))
    out.append("")

    lead_ph, lead_wh = groups(lead, PHARMACOPHORE), groups(lead, WARHEAD)
    if ref:
        rm = Chem.MolFromSmiles(ref["smiles"])
        if rm is not None:
            ref_ph, ref_wh = groups(rm, PHARMACOPHORE), groups(rm, WARHEAD)
            mcs = shared_core(lead, rm)
            out.append("PUBLISHED COMPOUND FOR THIS TARGET")
            out.append("  %s" % ref["name"][:80])
            out.append("  %s" % ref["smiles"])
            out.append("  %d heavy atoms, MW %.0f"
                       % (rm.GetNumHeavyAtoms(), Descriptors.MolWt(rm)))
            out.append("")
            out.append("SHARED CORE")
            if mcs["atoms"]:
                out.append("  %d atoms common to both (%.0f%% of the lead, %.0f%% of the published "
                           "compound)" % (mcs["atoms"],
                                          100.0 * mcs["atoms"] / lead.GetNumHeavyAtoms(),
                                          100.0 * mcs["atoms"] / rm.GetNumHeavyAtoms()))
                out.append("  %s" % (mcs["smarts"] or "")[:84])
            else:
                out.append("  none: the two molecules share no complete ring system")
            out.append("")
            out.append("PHARMACOPHORE PARITY")
            for g in sorted(set(ref_ph) & set(lead_ph)):
                out.append("  KEPT     %s" % g)
            for g in sorted(set(lead_ph) - set(ref_ph)):
                out.append("  ADDED    %s" % g)
            for g in sorted(set(ref_ph) - set(lead_ph)):
                out.append("  MISSING  %s   (the published compound has this and the lead does not)"
                           % g)
            out.append("")
            out.append("COVALENT HANDLES")
            out.append("  lead      : %s" % (", ".join(lead_wh) if lead_wh else "none"))
            out.append("  published : %s" % (", ".join(ref_wh) if ref_wh else "none"))
            cov_rounds = [r for r in rows
                          if str(r.get("goal_feature") or "").lower().startswith("coval")]
            for r in cov_rounds:
                out.append("  round %s aimed a covalent bond at %s. No scoring function used in "
                           "this study forms a bond, so the margin reported for that round does "
                           "not measure the selectivity the warhead was chosen for."
                           % (r.get("round"), r.get("goal_residue")))
    else:
        out.append("no published compound on file for this case")
        out.append("")
        out.append("LEAD PHARMACOPHORES: %s" % (", ".join(lead_ph) or "none recognised"))
        out.append("LEAD WARHEADS      : %s" % (", ".join(lead_wh) or "none"))
    out.append("")
    return "\n".join(out)

def main() -> int:
    refs = best_reference()
    blocks = []
    for path in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(path).replace("panel_plan_", "").replace(".csv", "")
        rows = campaign(path)
        if not rows:
            continue
        blocks.append(report(case, rows, refs.get(case)))
    if not blocks:
        print("no campaigns on disk")
        return 1
    text = "\n".join(blocks)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)
    print("wrote %s" % os.path.relpath(OUT, HERE))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
