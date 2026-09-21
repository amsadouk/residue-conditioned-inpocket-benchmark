from __future__ import annotations

import json
import urllib.request
from typing import Dict, List, Optional, Tuple

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"

WARHEADS: List[Tuple[str, str, str]] = [
    ("acrylamide", "michael", "C=CC(=O)[NX3]"),
    ("propiolamide", "michael_yne", "C#CC(=O)[NX3]"),
    ("vinyl sulfone", "michael", "C=C[SX4](=O)(=O)"),
    ("chloroacetamide", "sn2", "[Cl]C[CX3](=O)[NX3]"),
    ("alpha-halo carbonyl", "sn2", "[F,Cl,Br]C[CX3]=O"),
    ("acrylonitrile", "michael_cn", "C=C[CX2]#[NX1]"),
    ("vinyl ketone", "michael", "C=C[CX3](=O)[#6]"),
    ("acrylate ester", "michael", "C=C[CX3](=O)[OX2]"),
]

_SEQ_CACHE: Dict[str, Optional[str]] = {}


def detect(smiles: str) -> List[Tuple[str, str]]:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return []
    out = []
    for name, mech, sma in WARHEADS:
        q = Chem.MolFromSmarts(sma)
        if q is not None and m.HasSubstructMatch(q):
            out.append((name, mech))
    return out


def _uniprot_sequence(acc: str) -> Optional[str]:
    if acc in _SEQ_CACHE:
        return _SEQ_CACHE[acc]
    seq = None
    try:
        raw = urllib.request.urlopen(UNIPROT_URL.format(acc=acc), timeout=60).read()
        seq = ((json.loads(raw.decode("utf-8")).get("sequence") or {}).get("value")) or None
    except Exception:
        seq = None
    _SEQ_CACHE[acc] = seq
    return seq


def _align_pair(seq_a: str, seq_b: str):
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        al = PairwiseAligner()
        al.substitution_matrix = substitution_matrices.load("BLOSUM62")
        al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
        best = al.align(seq_a, seq_b)[0]
        return str(best[0]), str(best[1])
    except Exception:
        return None


def _map_positions(seq_a: str, seq_b: str) -> Optional[Dict[int, Optional[int]]]:
    aln = _align_pair(seq_a, seq_b)
    if not aln:
        return None
    a_aln, b_aln = aln
    mapping: Dict[int, Optional[int]] = {}
    ia = ib = 0
    for ca, cb in zip(a_aln, b_aln):
        if ca != "-":
            ia += 1
        if cb != "-":
            ib += 1
        if ca != "-":
            mapping[ia] = ib if cb != "-" else None
    return mapping
