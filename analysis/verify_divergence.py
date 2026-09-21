from __future__ import annotations

import os
import sys
import urllib.request
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as A

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "_pdb_cache")

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

CASES: Dict[str, Dict] = {
    "ABL1_SRC":   {"pdb": "6NPV", "chain": "A", "target": "P00519", "anti": "P12931",
                   "anti_name": "SRC"},
    "CA2_CA1":    {"pdb": "6SBL", "chain": "A", "target": "P00918", "anti": "P00915",
                   "anti_name": "CA1"},
    "EGFR_T790M": {"pdb": "6TFV", "chain": "A", "target": "P00533", "anti": "P00533",
                   "anti_name": "EGFR wild type", "point_mutations": {790: ("T", "M")}},
    "KRAS_G12C":  {"pdb": "6OIM", "chain": "A", "target": "P01116", "anti": "P01116",
                   "anti_name": "KRAS wild type", "point_mutations": {12: ("G", "C")}},
    "AKT1_AKT2":  {"pdb": "7NH5", "chain": "A", "target": "P31749", "anti": "P31751",
                   "anti_name": "AKT2"},
    "PARP1_PARP2": {"pdb": "7AAC", "chain": "A", "target": "P09874", "anti": "Q9UGN5",
                    "anti_name": "PARP2"},
}

_THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
          "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
          "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
          "Y": "Tyr", "V": "Val"}

def three(one: str) -> str:
    return _THREE.get((one or "").upper(), one or "?")

def fetch(pdb_id: str) -> Optional[str]:
    if not os.path.isdir(CACHE):
        os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, "%s.pdb" % pdb_id.upper())
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return open(path, encoding="utf-8", errors="ignore").read()
    try:
        txt = urllib.request.urlopen(
            "https://files.rcsb.org/download/%s.pdb" % pdb_id.upper(), timeout=60
        ).read().decode("utf-8", "ignore")
    except Exception:
        return None
    open(path, "w", encoding="utf-8").write(txt)
    return txt

def observed_chain(pdb_text: str, chain: str) -> Dict[int, str]:
    out: Dict[int, str] = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) > 26 and line[21] == chain:
            name = line[17:20].strip()
            if name in AA3:
                try:
                    out[int(line[22:26])] = AA3[name]
                except ValueError:
                    continue
    return out

def check(case: str, resnum: int) -> Dict:
    cfg = CASES.get(case)
    if not cfg:
        return {"case": case, "resnum": resnum, "ok": False, "reason": "case not configured"}

    txt = fetch(cfg["pdb"])
    if not txt:
        return {"case": case, "resnum": resnum, "ok": False,
                "reason": "could not retrieve %s" % cfg["pdb"]}
    obs = observed_chain(txt, cfg["chain"])
    if resnum not in obs:
        return {"case": case, "resnum": resnum, "ok": False,
                "reason": "residue %d absent from %s chain %s" % (resnum, cfg["pdb"], cfg["chain"])}
    struct_res = obs[resnum]

    tseq = A._uniprot_sequence(cfg["target"])
    aseq = A._uniprot_sequence(cfg["anti"])
    if not tseq or not aseq:
        return {"case": case, "resnum": resnum, "ok": False, "reason": "UniProt fetch failed"}

    nums = sorted(obs)
    sseq = "".join(obs[n] for n in nums)
    s2u = A._map_positions(sseq, tseq) or {}
    upos = s2u.get(nums.index(resnum) + 1)
    if not upos:
        return {"case": case, "resnum": resnum, "ok": False,
                "reason": "could not map structure position to the UniProt entry"}
    target_res = tseq[upos - 1]

    muts = cfg.get("point_mutations") or {}
    if upos in muts:
        wt, mt = muts[upos]
        target_res = mt
        anti_res, anti_pos = wt, upos
    elif cfg["target"] == cfg["anti"]:
        anti_res, anti_pos = target_res, upos
    else:
        t2a = A._map_positions(tseq, aseq) or {}
        anti_pos = t2a.get(upos)
        anti_res = aseq[anti_pos - 1] if anti_pos else None

    divergent = bool(anti_res and anti_res != target_res)
    return {"case": case, "resnum": resnum, "ok": True, "pdb": cfg["pdb"],
            "structure_residue": struct_res, "uniprot_pos": upos, "target_residue": target_res,
            "anti_name": cfg["anti_name"], "anti_pos": anti_pos, "anti_residue": anti_res,
            "divergent": divergent,
            "statement": ("%s %s%d aligns to %s %s%s, which differs"
                          % (case.split("_")[0], target_res, upos, cfg["anti_name"],
                             anti_res or "?", anti_pos or "?")) if divergent else
                         ("%s %s%d is CONSERVED in %s"
                          % (case.split("_")[0], target_res, upos, cfg["anti_name"])),
            }

def main() -> int:
    warheads = [("ABL1_SRC", 388), ("CA2_CA1", 206), ("EGFR_T790M", 797), ("KRAS_G12C", 12)]
    print("COVALENT WARHEAD RESIDUES, CHECKED AGAINST THE ALIGNMENT")
    print()
    print("%-13s %-9s %-11s %-13s %-15s %s"
          % ("case", "structure", "canonical", "target", "anti-target", "verdict"))
    print("-" * 82)
    rows = []
    for case, num in warheads:
        r = check(case, num)
        rows.append(r)
        if not r.get("ok"):
            print("%-13s %-9s %s" % (case[:13], num, r.get("reason")))
            continue
        print("%-13s %-9s %-11s %-13s %-15s %s"
              % (case[:13], "%s %d" % (r["pdb"], num), "%s%d" % (r["target_residue"],
                                                                 r["uniprot_pos"]),
                 r["structure_residue"] + str(num),
                 "%s%s" % (r["anti_residue"], r["anti_pos"]),
                 "DIVERGENT" if r["divergent"] else "conserved"))
    print()
    div = [r for r in rows if r.get("divergent")]
    con = [r for r in rows if r.get("ok") and not r.get("divergent")]
    print("%d of %d warhead residues are divergent in the anti-target." % (len(div), len(rows)))
    for r in con:
        print("  %s is conserved: %s. Its selection is a reach and reactivity decision, not a "
              "sequence-divergence one, and must not be reported under the divergence claim."
              % (r["case"], r["statement"]))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
