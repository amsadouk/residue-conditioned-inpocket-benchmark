from __future__ import annotations

import csv
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import rdFMCS, Descriptors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from structure_vs_literature import PHARMACOPHORE, AUXILIARY, ROLES, is_anchor, role_note, groups
import common as CE

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
MASTER = os.path.join(os.path.dirname(HERE), "data", "tables", "master_results.csv")
OUT = os.path.join(HERE, "TABLE_vs_commercial")

NAMED = re.compile(
    r"(imatinib|nilotinib|dasatinib|ponatinib|bosutinib|staurosporine|vemurafenib|dabrafenib|"
    r"encorafenib|regorafenib|sorafenib|sb590885|palbociclib|ribociclib|abemaciclib|gefitinib|"
    r"erlotinib|afatinib|osimertinib|dacomitinib|neratinib|lapatinib|tofacitinib|ruxolitinib|"
    r"baricitinib|fedratinib|celecoxib|rofecoxib|indometacin|acetazolamide|methazolamide|"
    r"dorzolamide|brinzolamide|sotorasib|adagrasib|vandetanib|dup-697|olaparib|niraparib|"
    r"rucaparib|talazoparib|veliparib|alpelisib|idelalisib|duvelisib|copanlisib|crizotinib|"
    r"axitinib|sunitinib|pazopanib|nintedanib|midostaurin|gilteritinib)", re.I)

CASE_ALIAS = {
    "EGFR_T790M_panel": "EGFR_T790M",
    "EGFR_T790M_vs_WT": "EGFR_T790M",
}

ANTI_TARGET_RESIDUE = {
    ("ABL1_SRC", "Cys388"): "SRC does not carry a cysteine here",
    ("CA2_CA1", "Cys206"): "CA1 carries Ser206",
    ("EGFR_T790M", "Cys797"): "the reference also targets Cys797",
    ("KRAS_G12C", "Cys12"): "Cys12 is the G12C mutation itself",
}

APPROVED = {"imatinib", "tofacitinib", "vandetanib", "gefitinib", "erlotinib", "afatinib",
            "osimertinib", "dacomitinib", "neratinib", "lapatinib", "ruxolitinib", "baricitinib",
            "fedratinib", "celecoxib", "indometacin", "acetazolamide", "methazolamide",
            "dorzolamide", "brinzolamide", "sotorasib", "adagrasib", "vemurafenib", "dabrafenib",
            "encorafenib", "regorafenib", "sorafenib", "palbociclib", "ribociclib", "abemaciclib",
            "nilotinib", "dasatinib", "ponatinib", "bosutinib"}

def _f(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() in ("", "None", "nan", "-"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

MARKETED = os.path.join(os.path.dirname(HERE), "data", "tables", "marketed_drugs.csv")

def commercial() -> Dict[str, List[Dict]]:
    out: Dict[str, List[Dict]] = {}
    if os.path.exists(MARKETED):
        for r in csv.DictReader(open(MARKETED, newline="", encoding="utf-8")):
            if not r.get("smiles"):
                continue
            out.setdefault(r["case"], []).append(
                {"key": r["drug"], "name": r["drug"], "smiles": r["smiles"], "vina": None,
                 "approved": str(r.get("approved")).lower() in ("true", "1", "yes"),
                 "note": r.get("note") or ""})

    by_key = {(c, d["key"]): d for c, ds in out.items() for d in ds}
    for r in csv.DictReader(open(MASTER, newline="", encoding="utf-8")):
        if r.get("kind") not in ("reference", "drug") or not r.get("smiles"):
            continue
        nm = r.get("resolved_name") or r.get("name") or ""
        m = NAMED.search(nm)
        if not m:
            continue
        case = str(r.get("case") or "").split("__")[0]
        case = CASE_ALIAS.get(case, case)
        d = by_key.get((case, m.group(0).lower()))
        if d is None:
            continue
        a = _f(r.get("affinity_target"))
        if a is not None and (d["vina"] is None or a < d["vina"]):
            d["vina"] = a
    return out

def final_leads() -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    for f in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(f).replace("panel_plan_", "").replace(".csv", "")
        rows = [r for r in csv.DictReader(open(f, newline="", encoding="utf-8")) if r.get("smiles")]
        if not rows:
            continue
        rows.sort(key=lambda r: int(r.get("round") or 0))
        fin = rows[-1]
        cov = [r for r in rows if str(r.get("goal_feature") or "").lower().startswith("coval")]
        fin["_covalent_residues"] = [r.get("goal_residue") for r in cov]
        fin["_n_rounds"] = len(rows)
        out[case] = fin
    return out

def compare(lead: Dict, drug: Dict, case: str) -> Dict:
    L = Chem.MolFromSmiles(lead["smiles"])
    R = Chem.MolFromSmiles(drug["smiles"])
    if L is None or R is None:
        return {}
    try:
        mcs = rdFMCS.FindMCS([L, R], timeout=20, ringMatchesRingOnly=True,
                             completeRingsOnly=True)
        shared, pct = mcs.numAtoms, 100.0 * mcs.numAtoms / R.GetNumHeavyAtoms()
    except Exception:
        shared, pct = 0, 0.0
    lp, rp = set(groups(L, PHARMACOPHORE)), set(groups(R, PHARMACOPHORE))
    wh = [w[0] for w in CE.detect(lead["smiles"])]
    dwh = [w[0] for w in CE.detect(drug["smiles"])]
    res = lead.get("_covalent_residues") or []
    note = ""
    for r in res:
        note = ANTI_TARGET_RESIDUE.get((case, r), "")
        if note:
            break
    return {
        "case": case, "drug": drug["name"][:34], "drug_key": drug["key"],
        "approved": bool(drug.get("approved", drug["key"] in APPROVED)),
        "rounds": lead.get("_n_rounds"),
        "shared_atoms": shared, "shared_pct": round(pct, 0),
        "kept": sorted(g for g in (lp & rp) if is_anchor(g)),
        "added": sorted(g for g in (lp - rp) if is_anchor(g)),
        "missing": sorted(g for g in (rp - lp) if is_anchor(g)),
        "aux_added": sorted(g for g in (lp - rp) if not is_anchor(g)),
        "aux_dropped": sorted(g for g in (rp - lp) if not is_anchor(g)),
        "warhead": wh, "drug_warhead": dwh,
        "covalent_residue": (res[0] if res else None), "anti_target_note": note,
        "lead_vina": _f(lead.get("affinity_target")),
        "drug_vina": drug["vina"],
        "lead_ha": L.GetNumHeavyAtoms(), "drug_ha": R.GetNumHeavyAtoms(),
        "lead_smiles": lead["smiles"],
    }

def main() -> int:
    drugs, leads = commercial(), final_leads()
    rows: List[Dict] = []
    for case, lead in sorted(leads.items()):
        for d in drugs.get(case, []):
            c = compare(lead, d, case)
            if c:
                rows.append(c)
    if not rows:
        print("no case pairs a final lead with a named commercial compound")
        return 1

    print("FINAL-ROUND LEAD vs NAMED COMMERCIAL COMPOUND")
    print("structure decides; the docking score is context in the last column")
    print()
    hdr = ("%-15s %-15s %-4s %-9s %-4s %-24s %s"
           % ("case", "compound", "appr", "shared", "miss", "covalent warhead", "Vina lead/drug"))
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        wh = (r["warhead"][0] if r["warhead"] else "-")
        if r["warhead"] and r["covalent_residue"]:
            wh = "%s at %s" % (r["warhead"][0][:14], r["covalent_residue"])
        print("%-15s %-15s %-4s %-9s %-4s %-24s %s / %s"
              % (r["case"][:15], r["drug_key"][:15], "yes" if r["approved"] else "-",
                 "%d%%" % r["shared_pct"], len(r["missing"]), wh[:24],
                 r["lead_vina"], r["drug_vina"]))

    print()
    print("RECOGNITION CHEMISTRY, per pair")
    for r in rows:
        print("  %s vs %s" % (r["case"], r["drug_key"]))
        print("     kept    : %s" % (", ".join(r["kept"]) or "none"))
        print("     added   : %s" % (", ".join(r["added"]) or "none"))
        print("     MISSING : %s" % (", ".join(r["missing"]) or "NONE"))
        for g in r["aux_dropped"]:
            print("     not carried: %s [%s]" % (g, role_note(g) or "peripheral"))
        for g in r["aux_added"]:
            print("     added, peripheral: %s [%s]" % (g, role_note(g) or "peripheral"))
        if r["warhead"]:
            print("     COVALENT: %s at %s. %s. The compound has %s."
                  % (", ".join(r["warhead"]), r["covalent_residue"] or "a divergent cysteine",
                     r["anti_target_note"] or "the anti-target differs at this position",
                     (", ".join(r["drug_warhead"]) if r["drug_warhead"] else "no warhead")))

    n_nomiss = sum(1 for r in rows if not r["missing"])
    n_appr = sum(1 for r in rows if r["approved"])
    n_wh = sum(1 for r in rows if r["warhead"] and not r["drug_warhead"])
    print()
    print("SUMMARY")
    print("  pairs compared                                    %d" % len(rows))
    print("  against approved drugs                            %d" % n_appr)
    print("  leads losing NO recognition chemistry             %d" % n_nomiss)
    print("  leads adding a warhead the compound lacks         %d" % n_wh)
    print("  median shared core, as %% of the compound          %d%%"
          % sorted(r["shared_pct"] for r in rows)[len(rows) // 2])

    cols = ["case", "drug", "drug_key", "approved", "rounds", "shared_atoms", "shared_pct",
            "kept", "added", "missing", "aux_added", "aux_dropped", "warhead", "drug_warhead",
            "covalent_residue", "anti_target_note", "lead_vina", "drug_vina",
            "lead_ha", "drug_ha", "lead_smiles"]
    with open(OUT + ".csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (", ".join(v) if isinstance(v, list) else v) for k, v in r.items()})
    print("\nwrote %s.csv" % os.path.relpath(OUT, HERE))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
