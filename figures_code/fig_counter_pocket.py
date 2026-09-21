from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import verify_divergence as VD
import common as A

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "FIG_counter_pocket")

PANELS = [("ABL1_SRC", 388), ("CA2_CA1", 206), ("KRAS_G12C", 12), ("EGFR_T790M", 797)]

C_DIV = (0.20, 0.60, 0.33)
C_CON = (0.72, 0.74, 0.76)
C_WAR = (0.95, 0.55, 0.15)
C_RIB = (0.88, 0.89, 0.90)

POCKET_R = 11.0

def residue_atoms(pdb_text: str, chain: str) -> Dict[int, List[Tuple[str, np.ndarray]]]:
    out: Dict[int, List[Tuple[str, np.ndarray]]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or len(line) < 55 or line[21] != chain:
            continue
        if line[16] not in (" ", "A"):
            continue
        try:
            num = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        out.setdefault(num, []).append((line[12:16].strip(), xyz))
    return out

def ligand_centre(pdb_text: str, chain: str) -> Optional[np.ndarray]:
    groups: Dict[str, List[np.ndarray]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM") or len(line) < 55:
            continue
        resn = line[17:20].strip()
        if resn in ("HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "ZN", "MG", "NA", "CL", "CA"):
            continue
        if line[21] != chain:
            continue
        try:
            groups.setdefault("%s%s" % (resn, line[22:26]), []).append(
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        except ValueError:
            continue
    if not groups:
        return None
    best = max(groups.values(), key=len)
    return np.mean(np.array(best), axis=0)

def build(case: str, warhead_num: int) -> Optional[Dict]:
    cfg = VD.CASES.get(case)
    if not cfg:
        return None
    txt = VD.fetch(cfg["pdb"])
    if not txt:
        return None
    chain = cfg["chain"]
    atoms = residue_atoms(txt, chain)
    if not atoms:
        return None

    centre = ligand_centre(txt, chain)
    if centre is None:
        wa = atoms.get(warhead_num)
        if not wa:
            return None
        centre = np.mean(np.array([x for _, x in wa]), axis=0)

    pocket = []
    for num, ats in atoms.items():
        coords = np.array([x for _, x in ats])
        if np.min(np.linalg.norm(coords - centre, axis=1)) <= POCKET_R:
            side = [x for n, x in ats if n not in ("N", "C", "O")]
            pocket.append((num, np.mean(np.array(side or coords), axis=0)))

    tseq = A._uniprot_sequence(cfg["target"])
    aseq = A._uniprot_sequence(cfg["anti"])
    if not tseq or not aseq:
        return None
    obs = VD.observed_chain(txt, chain)
    obs_nums = sorted(obs)
    sseq = "".join(obs[k] for k in obs_nums)
    s2u = A._map_positions(sseq, tseq) or {}
    t2a = (A._map_positions(tseq, aseq) or {}) if cfg["target"] != cfg["anti"] else None
    muts = cfg.get("point_mutations") or {}

    def divergence(num: int):
        if num not in obs:
            return None
        upos = s2u.get(obs_nums.index(num) + 1)
        if not upos:
            return None
        tres = tseq[upos - 1]
        if upos in muts:
            wt, mt = muts[upos]
            return {"uniprot_pos": upos, "target_residue": mt, "anti_residue": wt,
                    "divergent": True}
        if t2a is None:
            return {"uniprot_pos": upos, "target_residue": tres, "anti_residue": tres,
                    "divergent": False}
        apos = t2a.get(upos)
        ares = aseq[apos - 1] if apos else None
        return {"uniprot_pos": upos, "target_residue": tres, "anti_residue": ares,
                "divergent": bool(ares and ares != tres)}

    rows = []
    for num, xyz in sorted(pocket):
        chk = divergence(num)
        if not chk:
            continue
        rows.append({"num": num, "xyz": xyz, "divergent": bool(chk.get("divergent")),
                     "target_residue": chk.get("target_residue"),
                     "anti_residue": chk.get("anti_residue"),
                     "uniprot_pos": chk.get("uniprot_pos"),
                     "is_warhead": (num == warhead_num)})
    ca = np.array([x for n, ats in sorted(atoms.items()) for an, x in ats if an == "CA"])
    return {"case": case, "pdb": cfg["pdb"], "anti_name": cfg["anti_name"],
            "rows": rows, "ca": ca, "centre": centre}

def render(d: Dict, path: str, size=(1100, 950)) -> Optional[str]:
    import pyvista as pv
    pv.OFF_SCREEN = True
    pl = pv.Plotter(off_screen=True, window_size=list(size))
    pl.set_background("white")

    ca = d["ca"]
    near = ca[np.linalg.norm(ca - d["centre"], axis=1) <= 26.0]
    if len(near) > 3:
        try:
            spl = pv.Spline(near, max(len(near) * 3, 60)).tube(radius=0.28)
            pl.add_mesh(spl, color=C_RIB, smooth_shading=True)
        except Exception:
            pass

    for r in d["rows"]:
        col = C_WAR if r["is_warhead"] else (C_DIV if r["divergent"] else C_CON)
        rad = 1.30 if (r["is_warhead"] or r["divergent"]) else 0.95
        pl.add_mesh(pv.Sphere(radius=rad, center=r["xyz"]), color=col, smooth_shading=True)
        if r["divergent"] or r["is_warhead"]:
            lab = "%s%d (UniProt %d)" % (r["target_residue"], r["num"], r["uniprot_pos"])
            if r["divergent"]:
                lab += " / %s" % r["anti_residue"]
            pl.add_point_labels([r["xyz"]], [lab], font_size=17, text_color="black",
                                shape=None, always_visible=True, show_points=False)

    pl.camera_position = "xy"
    pl.camera.zoom(1.25)
    try:
        pl.screenshot(path)
    except Exception as exc:
        print("  render failed: %s" % str(exc)[:90])
        return None
    finally:
        pl.close()
    return path if os.path.exists(path) else None

def main() -> int:
    S.use()
    shots = []
    for case, num in PANELS:
        d = build(case, num)
        if not d:
            print("  %s: could not build" % case)
            continue
        png = os.path.join(HERE, "_cp_%s.png" % case)
        if render(d, png):
            shots.append((d, png))
            nd = sum(1 for r in d["rows"] if r["divergent"])
            print("  %-14s %s  %d pocket residues, %d divergent, warhead %s"
                  % (case, d["pdb"], len(d["rows"]), nd,
                     "on a divergent position" if any(
                         r["is_warhead"] and r["divergent"] for r in d["rows"])
                     else "on a CONSERVED position"))
    if not shots:
        print("nothing rendered")
        return 1

    n = len(shots)
    cols = 2
    rowsn = (n + cols - 1) // cols
    fig = plt.figure(figsize=(S.COL_2, 3.05 * rowsn))
    gs = fig.add_gridspec(rowsn, cols, hspace=0.20, wspace=0.03,
                          left=0.01, right=0.99, top=0.945, bottom=0.055)

    for i, (d, png) in enumerate(shots):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        ax.imshow(mpimg.imread(png))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        nd = sum(1 for r in d["rows"] if r["divergent"])
        wh = next((r for r in d["rows"] if r["is_warhead"]), None)
        ax.set_title(S.clean("%s  against %s" % (d["case"].split("_")[0], d["anti_name"])),
                     fontsize=7.2, loc="left", pad=3, color=S.C_TEXT, fontweight="bold")
        sub = "%d of %d pocket positions differ" % (nd, len(d["rows"]))
        if wh is not None:
            sub += ".  warhead at %s%d, %s" % (
                wh["target_residue"], wh["uniprot_pos"],
                "which differs" if wh["divergent"] else "which is CONSERVED")
        ax.text(0.0, -0.035, S.clean(sub), transform=ax.transAxes, fontsize=5.5,
                va="top", ha="left",
                color=("#2a7d4f" if (wh is not None and wh["divergent"]) else "#8c8c8c"))

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), dpi=450 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf and .png  (%d panels)" % (os.path.relpath(OUT, HERE), len(shots)))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
