from __future__ import annotations

import csv
import glob
import itertools
import os
import random
import sys
from typing import Dict, List, Optional, Tuple

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from rdkit.Chem import rdFMCS

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_structural import mcs_pct, load, decoys, N_DECOY, SEED

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "STATS_chemotype_scale.txt")

def best_against(mol, drug_mols: List) -> Tuple[float, Optional[str]]:
    best, who = 0.0, None
    for name, d in drug_mols:
        v = mcs_pct(mol, d)
        if v > best:
            best, who = v, name
    return best, who

def main() -> int:
    rng = random.Random(SEED)
    leads, drugs, actives = load()
    lines: List[str] = []

    def say(s=""):
        print(s)
        lines.append(s)

    say("CHEMOTYPE RECOVERY ON THE SCALE THE MARKETED DRUGS THEMSELVES SET")
    say("No docking score is used anywhere in this file.")
    say()

    say("1  HOW WELL DO A TARGET'S OWN MARKETED DRUGS RECOVER EACH OTHER?")
    say("   Shared core between each pair of approved drugs for the same target. This is what")
    say("   'two molecules that both work on this target' looks like as a number.")
    say()
    say("   %-16s %-5s %-9s %-9s %s" % ("target", "pairs", "median", "range", "example"))
    say("   " + "-" * 74)
    scale: Dict[str, List[float]] = {}
    for case, ds in sorted(drugs.items()):
        mols = [(d["key"], Chem.MolFromSmiles(d["smiles"])) for d in ds]
        mols = [(n, m) for n, m in mols if m is not None]
        if len(mols) < 2:
            continue
        vals, ex = [], ""
        for (na, ma), (nb, mb) in itertools.combinations(mols, 2):
            v = mcs_pct(ma, mb)
            vals.append(v)
            if v == max(vals):
                ex = "%s/%s %.0f%%" % (na[:9], nb[:9], v)
        vals.sort()
        scale[case] = vals
        say("   %-16s %-5d %-9s %-9s %s"
            % (case[:16], len(vals), "%.0f%%" % vals[len(vals) // 2],
               "%.0f to %.0f%%" % (vals[0], vals[-1]), ex))
    allpairs = sorted(v for vs in scale.values() for v in vs)
    if allpairs:
        say()
        say("   ACROSS ALL TARGETS: %d drug pairs, median %.0f%%, 90th percentile %.0f%%."
            % (len(allpairs), allpairs[len(allpairs) // 2],
               allpairs[int(0.9 * (len(allpairs) - 1))]))
        say("   Two marketed drugs for the SAME target typically share about %.0f%% of one"
            % allpairs[len(allpairs) // 2])
        say("   another. That is the scale every number below should be read against.")

    say()
    say("2  DID THE LEAD LAND ON A MARKETED CHEMOTYPE?")
    say("   Statistic: the best shared core across that target's marketed drugs.")
    say("   Null: each decoy scored against the SAME drug set, its own best taken. Both sides")
    say("   get best-of-N, so the maximum is not an advantage handed to one side.")
    say()
    say("   %-16s %-8s %-14s %-11s %-13s %s"
        % ("target", "lead", "matched", "decoy null", "drug-to-drug", "p (exact)"))
    say("   " + "-" * 82)

    rows = []
    n_sig = 0
    for case, lead in sorted(leads.items()):
        L = Chem.MolFromSmiles(lead.get("smiles") or "")
        ds = drugs.get(case) or []
        if L is None or not ds:
            continue
        mols = [(d["key"], Chem.MolFromSmiles(d["smiles"])) for d in ds]
        mols = [(n, m) for n, m in mols if m is not None]
        if not mols:
            continue

        lead_best, who = best_against(L, mols)
        pool = [Chem.MolFromSmiles(s) for s in decoys(case, actives, rng)]
        pool = [m for m in pool if m is not None]
        null = [best_against(m, mols)[0] for m in pool]
        if not null:
            continue
        n_ge = sum(1 for x in null if x >= lead_best)
        p = (n_ge + 1) / float(len(null) + 1)
        med_null = sorted(null)[len(null) // 2]
        dd = scale.get(case) or []
        dd_med = (sorted(dd)[len(dd) // 2] if dd else None)
        if p < 0.05:
            n_sig += 1
        say("   %-16s %-8s %-14s %-11s %-13s %.4f%s"
            % (case[:16], "%.0f%%" % lead_best, (who or "")[:14],
               "%.0f%%" % med_null,
               ("%.0f%%" % dd_med) if dd_med is not None else "n/a",
               p, "  *" if p < 0.05 else ""))
        rows.append({"case": case, "lead_best": round(lead_best, 1), "matched_drug": who,
                     "null_median": round(med_null, 1),
                     "drug_to_drug_median": (round(dd_med, 1) if dd_med is not None else ""),
                     "p": round(p, 4),
                     "beats_drug_to_drug": bool(dd_med is not None and lead_best > dd_med)})

    say()
    say("   %d of %d targets beat their matched decoy null at p<0.05." % (n_sig, len(rows)))
    above = [r for r in rows if r["beats_drug_to_drug"]]
    say("   %d of %d leads recover their nearest marketed drug BETTER than two marketed drugs"
        % (len(above), len(rows)))
    say("   for that same target recover each other.")
    if rows:
        lb = sorted(r["lead_best"] for r in rows)
        say("   Median best recovery %.0f%%, against a drug-to-drug median of %.0f%%."
            % (lb[len(lb) // 2], allpairs[len(allpairs) // 2] if allpairs else float("nan")))

    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    if rows:
        with open(os.path.join(HERE, "STATS_chemotype_scale.csv"), "w", newline="",
                  encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            for r in rows:
                w.writerow(r)
    print("\nwrote STATS_chemotype_scale.txt and .csv")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
