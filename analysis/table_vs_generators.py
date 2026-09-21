from __future__ import annotations

import csv
import glob
import os
import sys
from typing import Dict, List

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
OUT = os.path.join(HERE, "TABLE_vs_generators")

GENERATORS: List[Dict] = [
    {"name": "liGAN", "year": 2022, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "voxel generative model; evaluated on 10 targets, not the 100-pocket split"},
    {"name": "AR (3D-SBDD)", "year": 2021, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "autoregressive atom placement, NeurIPS"},
    {"name": "GraphBP", "year": 2022, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "atom placement in local spherical frames, ICML"},
    {"name": "Pocket2Mol", "year": 2022, "ha": 17.7, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "reported mean size is below the smallest molecule generated here"},
    {"name": "FLAG", "year": 2023, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "motif-by-motif fragment assembly, ICLR"},
    {"name": "DrugGPS", "year": 2023, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "subpocket prototypes, generalises to unseen targets, ICML"},
    {"name": "SurfGen", "year": 2023, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "surface-mesh pocket encoder, Nat Comput Sci"},
    {"name": "TargetDiff", "year": 2023, "ha": 24.2, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "diffusion, pocket conditioned, ICLR"},
    {"name": "DecompDiff", "year": 2023, "ha": 29.4, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "decomposed priors; best redocked score at the largest size"},
    {"name": "ResGen", "year": 2023, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": None, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "Nat Mach Intell; full text inaccessible, selectivity UNVERIFIED"},
    {"name": "D3FG", "year": 2023, "ha": None, "pocket3d": True, "residue": "partial",
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "functional groups as rigid bodies; group-level, not residue-level, conditioning"},
    {"name": "Lingo3DMol", "year": 2023, "ha": None, "pocket3d": True, "residue": "partial",
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "fragment-SMILES LM that predicts non-covalent interaction sites"},
    {"name": "DiffSBDD", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": True, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "explicit negative design on ONE kinase pair, BIKE over MPSK1, as a case study"},
    {"name": "PMDM", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "dual-scale equivariant diffusion, Nat Commun"},
    {"name": "PocketFlow", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "chemical-knowledge constraints, wet-lab hits for HAT1 and YTHDC1, Nat Mach Intell"},
    {"name": "IPDiff", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "pretrained interaction prior shifts the diffusion trajectory, ICLR"},
    {"name": "IRDiff", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "retrieval-augmented diffusion steered by reference ligands, ICML"},
    {"name": "AliDiff", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "preference alignment on binding energy, NeurIPS"},
    {"name": "TacoGFN", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "target-conditioned GFlowNet, reward is affinity plus SA plus QED, TMLR"},
    {"name": "FlexSBDD", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "flow matching over a FLEXIBLE pocket, relaxes side chains, NeurIPS"},
    {"name": "DeepICL", "year": 2024, "ha": None, "pocket3d": True, "residue": True,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": "partial",
     "note": "THE residue-conditioning precedent: explicit interaction type per pocket atom"},
    {"name": "PIDiff", "year": 2024, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": "partial", "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "physics-informed; off-target checked AFTER generation, not optimised for"},
    {"name": "BindGPT", "year": 2025, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "one language model for 3D generation, RL fine-tuned, AAAI"},
    {"name": "MolPilot", "year": 2025, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "VLB-optimal scheduling; 95.9% PoseBusters pass, the geometry state of the art"},
    {"name": "DynamicFlow", "year": 2025, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "apo to holo full-atom flow, integrates protein dynamics, ICLR"},
    {"name": "CByG", "year": 2025, "ha": None, "pocket3d": True, "residue": False,
     "paralogue": True, "sel_in_search": True, "covalent": False, "rationale": False,
     "note": "CLOSEST ML NEIGHBOUR: selectivity as a guidance term plus an on/off-target "
             "benchmark; conditions on the pocket, not on a named residue"},
    {"name": "CMD-GEN", "year": 2025, "ha": None, "pocket3d": True, "residue": "partial",
     "paralogue": True, "sel_in_search": True, "covalent": False, "rationale": False,
     "note": "CLOSEST EXPERIMENTAL NEIGHBOUR: selective pharmacophore points; wet-lab "
             "PARP1/PARP2 selective compound at >787-fold. Stronger evidence than ours"},
    {"name": "DiffPharma", "year": 2026, "ha": None, "pocket3d": True, "residue": True,
     "paralogue": False, "sel_in_search": False, "covalent": False, "rationale": False,
     "note": "interaction-particle constraints at residue level, npj Drug Discovery"},
    {"name": "MCTS kinase sel.", "year": 2022, "ha": None, "pocket3d": False, "residue": False,
     "paralogue": True, "sel_in_search": True, "covalent": False, "rationale": False,
     "note": "multiobjective MCTS for selectivity among kinase homologs, JCIM 62:5351; "
             "the canonical selectivity prior art, but ligand-based"},
    {"name": "CogMol", "year": 2020, "ha": None, "pocket3d": False, "residue": False,
     "paralogue": True, "sel_in_search": True, "covalent": False, "rationale": False,
     "note": "off-target penalty in a SMILES-VAE latent space, not 3D"},
    {"name": "CovaGEN", "year": 2026, "ha": None, "pocket3d": False, "residue": False,
     "paralogue": False, "sel_in_search": False, "covalent": True, "rationale": False,
     "note": "THE covalent precedent: latent diffusion with warhead guidance, validated on "
             "EGFR T790M. Sequence-conditioned via ESM-2, so no pocket geometry"},
    {"name": "Covalent MORL", "year": 2026, "ha": None, "pocket3d": False, "residue": "partial",
     "paralogue": False, "sel_in_search": False, "covalent": True, "rationale": False,
     "note": "SMILES-LSTM policy gradient with a residue-affinity reward term"},
    {"name": "This work", "year": 2026, "ha": None, "pocket3d": True, "residue": True,
     "paralogue": True, "sel_in_search": True, "covalent": True, "rationale": True,
     "note": "same-frame counter-pocket enters multi-objective SELECTION; warhead placed at a "
             "named divergent cysteine; every round states its reason"},
]

def our_size() -> Dict[str, float]:
    ha = []
    for f in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        for r in csv.DictReader(open(f, newline="", encoding="utf-8")):
            m = Chem.MolFromSmiles(r.get("smiles") or "")
            if m is not None:
                ha.append(m.GetNumHeavyAtoms())
    if not ha:
        return {}
    ha.sort()
    return {"n": len(ha), "mean": sum(ha) / len(ha), "median": ha[len(ha) // 2],
            "min": ha[0], "max": ha[-1]}

def main() -> int:
    size = our_size()
    for g in GENERATORS:
        if g["name"] == "This work" and size:
            g["ha"] = round(size["mean"], 1)

    def mark(b):
        if b == "partial":
            return "partial"
        if b is None:
            return "unknown"
        return "yes" if b else "no"

    rows = []
    print("STRUCTURE-BASED GENERATORS, COMPARED ON CAPABILITY")
    print("Docking is not compared here; against these methods it is compared in "
          "FIG_sbdd_comparison.")
    print()
    hdr = ("%-17s %-5s %-8s %-9s %-10s %-9s %-9s %s"
           % ("method", "year", "3D", "residue", "paralogue", "sel in", "covalent", "stated"))
    print(hdr)
    print("%-17s %-5s %-8s %-9s %-10s %-9s %-9s %s"
          % ("", "", "pocket", "cond.", "anti-tgt", "search", "design", "reason"))
    print("-" * len(hdr))
    for g in GENERATORS:
        print("%-17s %-5s %-8s %-9s %-10s %-9s %-9s %s"
              % (g["name"][:17], g.get("year", ""), mark(g.get("pocket3d")),
                 mark(g["residue"]), mark(g["paralogue"]), mark(g["sel_in_search"]),
                 mark(g["covalent"]), mark(g["rationale"])))
        rows.append(g)

    print()
    others = [g for g in GENERATORS if g["name"] != "This work"]
    both = [g for g in others
            if g.get("pocket3d") and g["paralogue"] is True and g["covalent"] is True]
    p3_par = [g for g in others if g.get("pocket3d") and g["paralogue"] is True]
    p3_res = [g for g in others if g.get("pocket3d") and g["residue"] is True]
    p3_cov = [g for g in others if g.get("pocket3d") and g["covalent"] is True]
    print("Among %d published methods surveyed:" % len(others))
    print("  %d are 3D pocket conditioned AND evaluate a paralogue  (%s)"
          % (len(p3_par), ", ".join(g["name"] for g in p3_par) or "none"))
    print("  %d are 3D pocket conditioned AND condition on a named residue  (%s)"
          % (len(p3_res), ", ".join(g["name"] for g in p3_res) or "none"))
    print("  %d are 3D pocket conditioned AND place a covalent warhead  (%s)"
          % (len(p3_cov), ", ".join(g["name"] for g in p3_cov) or "none"))
    print("  %d do all three, which is the gap this work occupies  (%s)"
          % (len(both), ", ".join(g["name"] for g in both) or "none"))
    print()
    if size:
        print("Molecules generated in this study: n=%d, mean %.1f heavy atoms, median %d, "
              "range %d to %d." % (size["n"], size["mean"], size["median"],
                                   size["min"], size["max"]))
        print("They are larger than the diffusion models' reported means, and the docking score")
        print("rises with heavy-atom count. The docking comparison against these methods is in")
        print("FIG_sbdd_comparison and is reported with a size-binned panel for that reason.")
    print()

    with open(OUT + ".csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\nwrote %s.csv" % os.path.relpath(OUT, HERE))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
