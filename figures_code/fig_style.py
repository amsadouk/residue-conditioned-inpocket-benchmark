from __future__ import annotations

from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1.0 / 25.4
COL_1 = 89 * MM
COL_15 = 120 * MM
COL_2 = 183 * MM

C_GEN = "#1b6ca8"
C_REF = "#8c8c8c"
C_WIN = "#2a8a4a"
C_LOSE = "#b4472e"
C_RULE = "#d9d9d9"
C_TEXT = "#1a1a1a"

PALETTE = [C_GEN, C_REF, C_WIN, C_LOSE]

def use() -> None:
    try:
        import scienceplots
        plt.style.use(["science", "nature"])
    except Exception:
        plt.style.use("default")
        plt.rcParams.update({
            "font.family": "serif", "axes.spines.top": False, "axes.spines.right": False,
            "axes.grid": False, "xtick.direction": "in", "ytick.direction": "in",
        })
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 600,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.labelcolor": C_TEXT, "text.color": C_TEXT,
        "axes.edgecolor": "#444444", "axes.linewidth": 0.6,
        "xtick.color": "#444444", "ytick.color": "#444444",
        "legend.frameon": False,
        "axes.prop_cycle": plt.cycler(color=PALETTE),
        "text.usetex": False, "figure.facecolor": "white", "axes.facecolor": "white",
    })

def clean(s) -> str:
    if s is None:
        return ""
    out = str(s)
    for bad, good in (("—", "-"), ("–", "-"), ("−", "-"),
                      ("‘", "'"), ("’", "'"), ("“", '"'), ("”", '"')):
        out = out.replace(bad, good)
    return out

def figure(width: float = COL_1, height: float = None, **kw) -> Tuple:
    h = height if height is not None else width * 0.72
    fig, ax = plt.subplots(figsize=(width, h), constrained_layout=True, **kw)
    return fig, ax

def panel_label(ax, letter: str, dx: float = -0.16, dy: float = 1.04) -> None:
    ax.text(dx, dy, clean(letter), transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="left", color=C_TEXT)

def despine(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

def save(fig, path_no_ext: str) -> str:
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (path_no_ext, ext))
    plt.close(fig)
    return path_no_ext + ".pdf"
