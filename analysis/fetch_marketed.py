from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from typing import Dict, List, Optional

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "data", "tables", "marketed_drugs.csv")
CACHE = os.path.join(HERE, "_pubchem_cache.json")

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/IsomericSMILES,InChIKey/JSON"

DRUGS: Dict[str, List[Dict]] = {
    "ABL1_SRC": [{"name": "imatinib", "approved": True}, {"name": "nilotinib", "approved": True},
                 {"name": "dasatinib", "approved": True}, {"name": "ponatinib", "approved": True},
                 {"name": "bosutinib", "approved": True}],
    "BRAF_V600E": [{"name": "vemurafenib", "approved": True},
                   {"name": "dabrafenib", "approved": True},
                   {"name": "encorafenib", "approved": True}],
    "BRAF_V600E_alt": [{"name": "vemurafenib", "approved": True},
                       {"name": "dabrafenib", "approved": True},
                       {"name": "encorafenib", "approved": True}],
    "CA2_CA1": [{"name": "acetazolamide", "approved": True},
                {"name": "methazolamide", "approved": True},
                {"name": "dorzolamide", "approved": True},
                {"name": "brinzolamide", "approved": True}],
    "CDK6_CDK2": [{"name": "palbociclib", "approved": True},
                  {"name": "ribociclib", "approved": True},
                  {"name": "abemaciclib", "approved": True}],
    "EGFR_T790M": [{"name": "gefitinib", "approved": True},
                   {"name": "erlotinib", "approved": True},
                   {"name": "afatinib", "approved": True},
                   {"name": "osimertinib", "approved": True},
                   {"name": "dacomitinib", "approved": True},
                   {"name": "neratinib", "approved": True}],
    "JAK2_JAK1": [{"name": "ruxolitinib", "approved": True},
                  {"name": "baricitinib", "approved": True},
                  {"name": "fedratinib", "approved": True},
                  {"name": "tofacitinib", "approved": True}],
    "KRAS_G12C": [{"name": "sotorasib", "approved": True},
                  {"name": "adagrasib", "approved": True}],
    "PTGS2_PTGS1": [{"name": "celecoxib", "approved": True},
                    {"name": "rofecoxib", "approved": False,
                     "note": "withdrawn from market in 2004"},
                    {"name": "etoricoxib", "approved": True},
                    {"name": "valdecoxib", "approved": False,
                     "note": "withdrawn from market in 2005"}],
    "PIK3CA_PIK3CD": [{"name": "alpelisib", "approved": True},
                      {"name": "idelalisib", "approved": True},
                      {"name": "duvelisib", "approved": True},
                      {"name": "copanlisib", "approved": True}],
    "PARP1_PARP2": [{"name": "olaparib", "approved": True},
                    {"name": "niraparib", "approved": True},
                    {"name": "rucaparib", "approved": True},
                    {"name": "talazoparib", "approved": True}],
    "AKT1_AKT2": [{"name": "capivasertib", "approved": True},
                  {"name": "ipatasertib", "approved": False,
                   "note": "investigational, phase 3"}],
}

def _cache() -> Dict:
    if os.path.exists(CACHE):
        try:
            return json.load(open(CACHE, encoding="utf-8"))
        except Exception:
            return {}
    return {}

def resolve(name: str, cache: Dict) -> Optional[Dict]:
    if name in cache:
        return cache[name]
    try:
        url = PUG % urllib.parse.quote(name)
        raw = urllib.request.urlopen(url, timeout=45).read().decode("utf-8", "ignore")
        props = (json.loads(raw).get("PropertyTable") or {}).get("Properties") or []
        if not props:
            return None
        p = props[0]
        smi = (p.get("SMILES") or p.get("IsomericSMILES")
               or p.get("ConnectivitySMILES") or p.get("CanonicalSMILES"))
        if not smi or Chem.MolFromSmiles(smi) is None:
            return None
        rec = {"smiles": smi, "inchikey": p.get("InChIKey"), "cid": p.get("CID")}
        cache[name] = rec
        time.sleep(0.25)
        return rec
    except Exception:
        return None

def main() -> int:
    cache = _cache()
    rows, failed = [], []
    for case, drugs in sorted(DRUGS.items()):
        for d in drugs:
            rec = resolve(d["name"], cache)
            if not rec:
                failed.append("%s/%s" % (case, d["name"]))
                continue
            m = Chem.MolFromSmiles(rec["smiles"])
            rows.append({"case": case, "drug": d["name"], "approved": d.get("approved", False),
                         "note": d.get("note", ""), "smiles": rec["smiles"],
                         "inchikey": rec.get("inchikey"), "cid": rec.get("cid"),
                         "heavy_atoms": m.GetNumHeavyAtoms()})
    json.dump(cache, open(CACHE, "w", encoding="utf-8"), indent=1)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    by_case: Dict[str, int] = {}
    for r in rows:
        by_case[r["case"]] = by_case.get(r["case"], 0) + 1
    print("MARKETED COMPOUNDS RESOLVED FROM PUBCHEM BY NAME")
    print()
    for c in sorted(by_case):
        names = ", ".join(r["drug"] for r in rows if r["case"] == c)
        print("  %-16s %d   %s" % (c, by_case[c], names))
    print()
    print("%d compounds across %d receptor pairs -> %s"
          % (len(rows), len(by_case), os.path.basename(OUT)))
    if failed:
        print("could not resolve (dropped, not guessed): %s" % ", ".join(failed))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
