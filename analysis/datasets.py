from __future__ import annotations

import csv
import glob
import json
import os
import sys
from typing import Dict, List

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RERUN = os.path.join(REPO, "data", "campaigns")

sys.path.insert(0, HERE)

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def campaigns() -> Dict:
    files = sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv")))
    with_leads, rounds, empty = [], {}, []
    for f in files:
        case = os.path.basename(f).replace("panel_plan_", "").replace(".csv", "")
        rows = [r for r in csv.DictReader(open(f, newline="", encoding="utf-8"))
                if (r.get("smiles") or "").strip()]
        if rows:
            with_leads.append(case)
            rounds[case] = len(rows)
        else:
            empty.append(case)
    return {"attempted": len(files), "with_leads": with_leads, "rounds": rounds, "empty": empty}

def marketed() -> Dict:
    p = os.path.join(HERE, "marketed_drugs.csv")
    if not os.path.exists(p):
        return {}
    rows = list(csv.DictReader(open(p, newline="", encoding="utf-8")))
    by_case: Dict[str, int] = {}
    approved = 0
    for r in rows:
        by_case[r["case"]] = by_case.get(r["case"], 0) + 1
        if str(r.get("approved")).lower() in ("true", "1", "yes"):
            approved += 1
    return {"n": len(rows), "cases": by_case, "approved": approved}

def sbdd_population() -> Dict:
    try:
        from fig_sbdd_comparison import octara
        rows = octara()
    except Exception as exc:
        return {"error": str(exc)[:90]}
    ha = sorted(r["ha"] for r in rows)
    return {"n": len(rows), "median_ha": ha[len(ha) // 2] if ha else None,
            "min_ha": ha[0] if ha else None, "max_ha": ha[-1] if ha else None,
            "with_strain": sum(1 for r in rows if r.get("strain") is not None),
            "with_pose_valid": sum(1 for r in rows if r.get("pv") in ("true", "false"))}

def covalent() -> Dict:
    from structure_vs_literature import WARHEAD
    goals, kept = [], []
    for f in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(f).replace("panel_plan_", "").replace(".csv", "")
        rows = [r for r in csv.DictReader(open(f, newline="", encoding="utf-8"))
                if (r.get("smiles") or "").strip()]
        if not rows:
            continue
        rows.sort(key=lambda r: int(r.get("round") or 0))
        if any(str(r.get("goal_feature") or "").lower().startswith("coval") for r in rows):
            goals.append(case)
        m = Chem.MolFromSmiles(rows[-1]["smiles"])
        if m is None:
            continue
        for name, sma in WARHEAD:
            q = Chem.MolFromSmarts(sma)
            if q is not None and m.HasSubstructMatch(q):
                kept.append((case, name.split(" (")[0]))
                break
    return {"with_covalent_goal": goals, "final_lead_carries_warhead": kept}

def gui() -> Dict:
    p = os.path.join(HERE, "gui_rounds", "manifest.json")
    if not os.path.exists(p):
        return {}
    m = json.load(open(p, encoding="utf-8"))
    return {"case": "CA2_CA1", "rounds": [r["round"] for r in m],
            "heavy_atoms": [r["heavy_atoms"] for r in m]}

def main() -> int:
    C = campaigns()
    M = marketed()
    B = sbdd_population()
    V = covalent()
    G = gui()

    lines: List[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("POPULATIONS USED IN THIS STUDY")
    say("Every count below is recomputed from the files the figures read.")
    say()

    say("1  CAMPAIGNS")
    say("   attempted                     %d receptor pairs" % C["attempted"])
    say("   produced a final lead         %d" % len(C["with_leads"]))
    if C["empty"]:
        say("   produced none                 %d (%s)" % (len(C["empty"]),
                                                          ", ".join(C["empty"])))
    say("   rounds per campaign           %s"
        % ", ".join("%s %d" % (c, C["rounds"][c]) for c in sorted(C["rounds"])))
    say("   USED BY: containment, trajectories, structure panel")
    say()

    say("2  MARKETED DRUG SET")
    if M:
        say("   compounds                     %d, of which %d approved" % (M["n"], M["approved"]))
        say("   receptor pairs covered        %d of %d campaigns with leads"
            % (len(M["cases"]), len(C["with_leads"])))
        missing = sorted(set(C["with_leads"]) - set(M["cases"]))
        if missing:
            say("   campaigns with NO marketed drug for their target: %s" % ", ".join(missing))
        say("   per pair                      %s"
            % ", ".join("%s %d" % (k, v) for k, v in sorted(M["cases"].items())))
    say("   USED BY: commercial comparison, chemotype scale, recovery and retention statistics")
    say()

    say("3  SBDD BENCHMARK POPULATION")
    if "error" in B:
        say("   unavailable: %s" % B["error"])
    else:
        say("   molecules                     %d, screened and de-duplicated" % B["n"])
        say("   heavy atoms                   median %s, range %s to %s"
            % (B["median_ha"], B["min_ha"], B["max_ha"]))
        say("   with a strain value           %d" % B["with_strain"])
        say("   with a pose-validity verdict  %d" % B["with_pose_valid"])
    say("   NOT the same population as 1: this is every generated molecule that passed the")
    say("   chemistry screen, not one final lead per campaign.")
    say("   USED BY: the generative-model comparison figure")
    say()

    say("4  COVALENT SUBSETS, WHICH ARE TWO DIFFERENT THINGS")
    say("   campaigns setting a covalent goal      %d (%s)"
        % (len(V["with_covalent_goal"]), ", ".join(V["with_covalent_goal"])))
    say("   final leads actually carrying a warhead %d (%s)"
        % (len(V["final_lead_carries_warhead"]),
           ", ".join("%s: %s" % (c, w) for c, w in V["final_lead_carries_warhead"])))
    say("   Reporting the first number as though it were the second credits campaigns with")
    say("   warheads their final molecules do not have.")
    say()

    say("5  GUI IN-POCKET VIEWS")
    if G:
        say("   one campaign, %s, rounds %s at %s heavy atoms"
            % (G["case"], G["rounds"], G["heavy_atoms"]))
        say("   These poses were REDOCKED for the figure; they are not the poses scored during")
        say("   the campaign, and the docking values are not shown.")
    say()

    out = os.path.join(HERE, "DATASETS.txt")
    open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("wrote %s" % os.path.basename(out))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
