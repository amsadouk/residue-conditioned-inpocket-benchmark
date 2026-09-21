from __future__ import annotations

import csv
import glob
import os
import subprocess
import sys
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as S

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(os.path.dirname(HERE), "data", "campaigns")
TABLE = os.path.join(os.path.dirname(HERE), "data", "tables", "TABLE_vs_commercial.csv")
OUT = os.path.join(HERE, "FIG_benchmark_table")

C_BAR = "#1b6ca8"
C_WIN = "#2a7d4f"
C_WARN = "#b4472e"
C_GREY = "#8c8c8c"

def refresh() -> None:
    try:
        subprocess.run([sys.executable, os.path.join(HERE, "table_vs_commercial.py")],
                       cwd=HERE, capture_output=True, timeout=900)
    except Exception:
        pass

def rows() -> List[Dict]:
    if not os.path.exists(TABLE):
        return []
    out = []
    for r in csv.DictReader(open(TABLE, newline="", encoding="utf-8")):
        try:
            r["shared_pct"] = float(r.get("shared_pct") or 0)
        except ValueError:
            r["shared_pct"] = 0.0
        r["missing_list"] = [x for x in (r.get("missing") or "").split(", ") if x]
        r["kept_list"] = [x for x in (r.get("kept") or "").split(", ") if x]
        r["approved"] = str(r.get("approved")).lower() in ("true", "yes", "1")
        out.append(r)
    out.sort(key=lambda r: (-r["shared_pct"], r["case"]))
    return out

def panel_state() -> Dict[str, int]:
    out = {}
    for f in sorted(glob.glob(os.path.join(RERUN, "panel_plan_*.csv"))):
        case = os.path.basename(f).replace("panel_plan_", "").replace(".csv", "")
        with open(f, newline="", encoding="utf-8") as fh:
            out[case] = max(0, sum(1 for _ in fh) - 1)
    return out

def main() -> int:
    refresh()
    R = rows()
    if not R:
        print("no comparisons available yet")
        return 1
    S.use()

    n = len(R)
    fig = plt.figure(figsize=(S.COL_2, 0.30 * n + 1.65))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(-2.6, n + 2.0)
    ax.axis("off")

    X = {"case": 1.0, "drug": 14.0, "bar": 28.5, "pct": 42.0,
         "kept": 45.0, "miss": 61.0, "cov": 74.5, "vina": 99.0}
    BARW = 12.0

    ax.text(X["case"], n + 1.35, S.clean("Every receptor against its named commercial compounds"),
            fontsize=8.0, fontweight="bold", color=S.C_TEXT)
    ax.text(X["case"], n + 0.85,
            S.clean("final carried round of each campaign; structure decides, the docking score "
                    "is context"),
            fontsize=5.6, color=C_GREY)

    hdr = n + 0.20
    for k, lab in (("case", "receptor pair"), ("drug", "commercial compound"),
                   ("bar", "shared core"), ("kept", "recognition kept"),
                   ("miss", "recognition missing"), ("cov", "covalent warhead")):
        ax.text(X[k], hdr, S.clean(lab), fontsize=5.7, fontweight="bold", color=S.C_TEXT)
    ax.text(X["vina"], hdr, S.clean("Vina, context"), fontsize=5.7, fontweight="bold",
            color=C_GREY, ha="right")
    ax.plot([X["case"], 99.5], [hdr - 0.28, hdr - 0.28], linewidth=0.7, color="#444444")

    for i, r in enumerate(R):
        y = n - 1 - i
        ax.text(X["case"], y, S.clean(r["case"][:17]), fontsize=5.5, color=S.C_TEXT, va="center")
        nm = r["drug_key"][:14]
        ax.text(X["drug"], y, S.clean(nm), fontsize=5.5, va="center",
                color=(S.C_TEXT if r["approved"] else C_GREY),
                fontweight=("bold" if r["approved"] else "normal"))
        if r["approved"]:
            ax.text(X["drug"] + 9.2, y, S.clean("approved"), fontsize=4.4, va="center",
                    color=C_WIN)

        w = BARW * (r["shared_pct"] / 100.0)
        ax.add_patch(plt.Rectangle((X["bar"], y - 0.22), BARW, 0.44,
                                   facecolor="#eef0f2", edgecolor="none"))
        ax.add_patch(plt.Rectangle((X["bar"], y - 0.22), w, 0.44,
                                   facecolor=C_BAR, edgecolor="none"))
        ax.text(X["pct"], y, S.clean("%d%%" % r["shared_pct"]), fontsize=5.3, va="center",
                ha="right", color=S.C_TEXT)

        kept = ", ".join(x.split(" (")[0] for x in r["kept_list"][:2]) or "none"
        ax.text(X["kept"], y, S.clean(kept[:24]), fontsize=5.0, va="center", color="#444444")

        if r["missing_list"]:
            miss = ", ".join(x.split(" (")[0] for x in r["missing_list"][:2])
            ax.text(X["miss"], y, S.clean(miss[:20]), fontsize=5.0, va="center", color=C_WARN)
        else:
            ax.text(X["miss"], y, S.clean("none"), fontsize=5.0, va="center",
                    color=C_WIN, fontweight="bold")

        wh = (r.get("warhead") or "").split(", ")[0]
        if wh:
            same = bool((r.get("drug_warhead") or "").strip())
            res = r.get("covalent_residue") or ""
            wh_name = wh.split(" (")[0]
            txt = ("%s at %s" % (wh_name, res)) if res else ("%s, no residue" % wh_name)
            ax.text(X["cov"], y, S.clean(txt[:32]), fontsize=4.6, va="center",
                    color=(C_WIN if same else C_BAR), fontweight="bold")
            if same:
                ax.text(X["cov"] + 0.0, y - 0.34,
                        S.clean("same residue as the drug's warhead"), fontsize=4.1,
                        va="center", color=C_WIN)
        else:
            ax.text(X["cov"], y, S.clean("-"), fontsize=5.0, va="center", color=C_GREY)

        lv, dv = r.get("lead_vina") or "", r.get("drug_vina") or ""
        ax.text(X["vina"], y, S.clean("%s / %s" % (lv or "n/a", dv or "n/a")), fontsize=4.9,
                va="center", ha="right", color=C_GREY)

    n_appr = sum(1 for r in R if r["approved"])
    n_nomiss_rows = sum(1 for r in R if not r["missing_list"])
    n_nomiss = len({r["case"] for r in R if not r["missing_list"]})
    n_wh = len({r["case"] for r in R if (r.get("warhead") or "").strip()})
    n_cases = len({r["case"] for r in R})
    ax.plot([X["case"], 99.5], [-0.55, -0.55], linewidth=0.7, color="#444444")
    ax.text(X["case"], -1.05,
            S.clean("%d comparisons across %d receptor pairs, %d against approved drugs. "
                    "%d of the %d comparisons lose no recognition chemistry, spread over %d "
                    "receptor pairs. %d leads carry a covalent warhead."
                    % (len(R), n_cases, n_appr, n_nomiss_rows, len(R), n_nomiss, n_wh)),
            fontsize=5.6, color=S.C_TEXT, fontweight="bold")

    st = panel_state()
    done = [c for c, k in st.items() if k]
    ax.text(X["case"], -1.65,
            S.clean("panel state: %d receptor pairs carry molecules (%s). Cases still running "
                    "are excluded from this figure rather than shown partial."
                    % (len(done), ", ".join("%s %d" % (c.split('_')[0], k)
                                            for c, k in sorted(st.items()) if k))),
            fontsize=4.6, color=C_GREY)
    ax.text(X["case"], -2.20,
            S.clean("Recognition motifs are compared by role: solubilising and metabolic handles "
                    "are peripheral and are excluded from both the kept and missing columns."),
            fontsize=4.6, color=C_GREY)

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (OUT, ext), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote %s.pdf and .png  (%d comparisons)" % (os.path.relpath(OUT, HERE), len(R)))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
