from __future__ import annotations

import csv
import glob
import os
import random
import sys
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import rdFMCS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from structure_vs_literature import PHARMACOPHORE, groups, is_anchor
from table_vs_commercial import NAMED, CASE_ALIAS, APPROVED
import common as CE

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
MASTER = os.path.join(os.path.dirname(HERE), "data", "tables", "master_results.csv")
MARKETED = os.path.join(os.path.dirname(HERE), "data", "tables", "marketed_drugs.csv")
OUT = os.path.join(HERE, "STATS_structural.txt")

N_DECOY = 60
SEED = 20240517

def _f(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() in ("", "None", "nan", "-"):
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

def mcs_pct(a, b) -> float:
    try:
        r = rdFMCS.FindMCS([a, b], timeout=10, ringMatchesRingOnly=True, completeRingsOnly=True)
    except Exception:
        return 0.0
    if not r or r.numAtoms == 0:
        return 0.0
    return 100.0 * r.numAtoms / b.GetNumHeavyAtoms()

def load() -> Tuple[Dict[str, Dict], Dict[str, List[Dict]], Dict[str, List[str]]]:
    leads: Dict[str, Dict] = {}
    for f in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(f).replace("panel_plan_", "").replace(".csv", "")
        rows = [r for r in csv.DictReader(open(f, newline="", encoding="utf-8")) if r.get("smiles")]
        if not rows:
            continue
        rows.sort(key=lambda r: int(r.get("round") or 0))
        fin = rows[-1]
        fin["_covalent"] = [r.get("goal_residue") for r in rows
                            if str(r.get("goal_feature") or "").lower().startswith("coval")]
        fin["_scaffold"] = rows[0].get("scaffold_smiles")
        leads[case] = fin

    actives: Dict[str, List[str]] = {}
    for r in csv.DictReader(open(MASTER, newline="", encoding="utf-8")):
        if r.get("kind") not in ("reference", "drug") or not r.get("smiles"):
            continue
        case = CASE_ALIAS.get(str(r.get("case") or "").split("__")[0],
                              str(r.get("case") or "").split("__")[0])
        actives.setdefault(case, []).append(r["smiles"])

    drugs: Dict[str, List[Dict]] = {}
    if os.path.exists(MARKETED):
        for r in csv.DictReader(open(MARKETED, newline="", encoding="utf-8")):
            if not r.get("smiles"):
                continue
            drugs.setdefault(r["case"], []).append(
                {"key": r["drug"], "name": r["drug"], "smiles": r["smiles"],
                 "approved": str(r.get("approved")).lower() in ("true", "1", "yes"),
                 "note": r.get("note") or ""})
    return leads, drugs, actives

def decoys(case: str, actives: Dict[str, List[str]], rng) -> List[str]:
    pool = [s for c, v in actives.items() if c != case for s in v]
    rng.shuffle(pool)
    return pool[:N_DECOY]

def main() -> int:
    rng = random.Random(SEED)
    leads, drugs, actives = load()
    lines: List[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("STRUCTURAL STATISTICS. No docking score is used in any test below.")
    say()

    say("1  CHEMOTYPE RECOVERY AGAINST CHANCE")
    say("   null: the known actives of the OTHER targets, scored the same way against the same")
    say("   published compound. p is the exact fraction of decoys matching or beating the lead.")
    say()
    say("   %-15s %-14s %-8s %-22s %s" % ("case", "compound", "lead", "decoy null", "p (exact)"))
    say("   " + "-" * 76)

    recov, retain = [], []
    for case, lead in sorted(leads.items()):
        L = Chem.MolFromSmiles(lead["smiles"] or "")
        if L is None:
            continue
        dpool = [Chem.MolFromSmiles(s) for s in decoys(case, actives, rng)]
        dpool = [m for m in dpool if m is not None]
        for d in drugs.get(case, []):
            R = Chem.MolFromSmiles(d["smiles"] or "")
            if R is None or not dpool:
                continue
            lead_pct = mcs_pct(L, R)
            null = [mcs_pct(m, R) for m in dpool]
            n_ge = sum(1 for x in null if x >= lead_pct)
            p = (n_ge + 1) / float(len(null) + 1)
            med = sorted(null)[len(null) // 2]
            say("   %-15s %-14s %-8s %-22s %.4f%s"
                % (case[:15], d["key"][:14], "%.0f%%" % lead_pct,
                   "median %.0f%%, max %.0f%%" % (med, max(null)), p,
                   "  *" if p < 0.05 else ""))
            recov.append({"case": case, "drug": d["key"], "approved": d["approved"],
                          "lead_pct": lead_pct, "null_median": med, "p": p})

            rp = {g for g in groups(R, PHARMACOPHORE) if is_anchor(g)}
            lp = {g for g in groups(L, PHARMACOPHORE) if is_anchor(g)}
            holds = rp.issubset(lp)
            n_hold = sum(1 for m in dpool
                         if rp.issubset({g for g in groups(m, PHARMACOPHORE) if is_anchor(g)}))
            retain.append({"case": case, "drug": d["key"], "lead_holds": holds,
                           "decoy_rate": n_hold / float(len(dpool)), "n_motifs": len(rp)})

    say()
    sig = [r for r in recov if r["p"] < 0.05]
    say("   %d of %d comparisons beat the decoy null at p<0.05" % (len(sig), len(recov)))
    if recov:
        say("   median lead recovery %.0f%% against median null %.0f%%"
            % (sorted(r["lead_pct"] for r in recov)[len(recov) // 2],
               sorted(r["null_median"] for r in recov)[len(recov) // 2]))

    say()
    say("2  RECOGNITION RETENTION")
    say("   does the lead hold EVERY primary recognition motif of the compound, and how often")
    say("   does an unrelated active do the same")
    say()
    held = [r for r in retain if r["lead_holds"]]
    say("   leads holding all motifs: %d of %d" % (len(held), len(retain)))
    if retain:
        mean_rate = sum(r["decoy_rate"] for r in retain) / len(retain)
        say("   decoys holding all motifs: %.1f%% on average" % (100 * mean_rate))
        try:
            from scipy import stats as sps
            b = sps.binomtest(len(held), len(retain), max(mean_rate, 1e-6), alternative="greater")
            say("   binomial test against the decoy rate: p = %.5f%s"
                % (b.pvalue, "  *" if b.pvalue < 0.05 else ""))
        except Exception as exc:
            say("   (scipy unavailable: %s)" % str(exc)[:50])

    say()
    say("3  CONTAINMENT")
    ok = tot = 0
    for case, lead in sorted(leads.items()):
        L = Chem.MolFromSmiles(lead["smiles"] or "")
        S_ = Chem.MolFromSmiles(lead.get("_scaffold") or "")
        if L is None or S_ is None:
            continue
        tot += 1
        ok += bool(L.HasSubstructMatch(S_))
    say("   scaffold present whole in the final lead: %d of %d campaigns" % (ok, tot))
    if tot:
        import math
        z = 1.959963985
        ph = ok / float(tot)
        den = 1.0 + z * z / tot
        cen = (ph + z * z / (2 * tot)) / den
        half = z * math.sqrt(ph * (1 - ph) / tot + z * z / (4 * tot * tot)) / den
        say("   Wilson 95%% interval: %.0f%% to %.0f%% (n=%d)"
            % (100 * max(0.0, cen - half), 100 * min(1.0, cen + half), tot))

    say()
    say("4  COVALENT PLACEMENT AT DIVERGENT RESIDUES")
    import verify_divergence as VD

    cov = [(c, l) for c, l in sorted(leads.items()) if l.get("_covalent")]
    say("   campaigns with a covalent goal: %d" % len(cov))
    say()
    say("   %-15s %-20s %-9s %-11s %s"
        % ("case", "warhead on final lead", "residue", "anti-target", "divergent"))
    say("   " + "-" * 72)
    kept, divergent = [], []
    for case, lead in cov:
        wh = [w[0] for w in CE.detect(lead["smiles"] or "")]
        res = lead["_covalent"][0]
        has = bool(wh)
        try:
            num = int("".join(ch for ch in str(res) if ch.isdigit()))
        except ValueError:
            num = None
        chk = VD.check(case, num) if num is not None else {"ok": False}
        anti = ("%s%s" % (chk.get("anti_residue"), chk.get("anti_pos"))
                if chk.get("ok") else "not checked")
        verdict = ("DIVERGENT" if chk.get("divergent") else
                   "conserved" if chk.get("ok") else "unknown")
        say("   %-15s %-20s %-9s %-11s %s"
            % (case[:15], (wh[0].split(" (")[0] if has else "NONE")[:20], res, anti, verdict))
        if has:
            kept.append(case)
            if chk.get("divergent"):
                divergent.append(case)
    say()
    say("   %d of %d campaigns with a covalent goal carry a warhead on the FINAL lead."
        % (len(kept), len(cov)))
    say("   Of those %d, %d sit at a position the anti-target differs at." % (len(kept),
                                                                              len(divergent)))
    say("   A conserved position is not a failed design. Reach, geometry and reactivity are")
    say("   legitimate grounds for placing a warhead, and EGFR Cys797 is the residue the")
    say("   third-generation inhibitors bond. It is simply not evidence for the selectivity")
    say("   claim, so it is counted separately rather than folded in.")

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(HERE, "STATS_recovery.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recov[0]) if recov else ["case"])
        w.writeheader()
        for r in recov:
            w.writerow(r)
    print("\nwrote STATS_structural.txt and STATS_recovery.csv")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
