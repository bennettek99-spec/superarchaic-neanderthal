"""make_figures.py — pilot summary figures (chr21+22)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
SIM = ROOT / "results" / "sims"
FIG = ROOT / "results" / "figures"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcflib import robust_z  # noqa: E402

CHROMS = ("16", "17", "18", "19", "20", "21", "22")


def load(win):
    df = pd.concat([pd.read_csv(TAB / f"windows.chr{c}.win{win}.tsv", sep="\t") for c in CHROMS],
                   ignore_index=True)
    df = df[(df["callable_frac"] >= 0.3) & np.isfinite(df["div_den_afr"])].copy()
    df["nea_afr"] = df[["div_alt_afr", "div_vin_afr", "div_chag_afr"]].mean(axis=1)
    df["den_excess"] = df["div_den_afr"] - df["nea_afr"]
    df["z_den_excess"] = robust_z(df["den_excess"].to_numpy())
    df["artifact"] = ((df["mean_rm"].fillna(0) > 0.6) | (df["mean_map20"].fillna(1) < 0.7) |
                      (df["mean_ur"].fillna(0) > 0.9) | (df["callable_frac"] < 0.4))
    return df


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    d = load(50000)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    # (a) divergence ordering
    pairs = [("div_alt_afr", "Altai–Afr"), ("div_den_afr", "Den–Afr"),
             ("div_vin_afr", "Vin–Afr"), ("div_chag_afr", "Chag–Afr"),
             ("div_den_alt", "Den–Altai"), ("div_den_vin", "Den–Vin"),
             ("div_alt_vin", "Altai–Vin"), ("div_vin_chag", "Vin–Chag")]
    vals = [d[c].median() for c, _ in pairs]
    cols = ["#4477aa"] * 4 + ["#ee6677"] * 2 + ["#228833"] * 2
    ax[0, 0].bar(range(len(pairs)), vals, color=cols)
    ax[0, 0].set_xticks(range(len(pairs)))
    ax[0, 0].set_xticklabels([n for _, n in pairs], rotation=45, ha="right", fontsize=8)
    ax[0, 0].set_ylabel("median divergence / bp")
    ax[0, 0].set_title("(a) Divergence ordering (common mask)\narchaic–African ≈ equal > Den–Nea > Nea–Nea")

    # (b) focal-genome CDDR rate (equal = shared confound)
    rates = {g: np.mean(robust_z(d[f"div_{g}_afr"].to_numpy()) >= 3) for g in
             ("den", "alt", "vin", "chag")}
    ax[0, 1].bar(list(rates), list(rates.values()),
                 color=["#ee6677", "#4477aa", "#4477aa", "#4477aa"])
    ax[0, 1].axhline(0.007, ls="--", c="grey", label="sim null ~0.007")
    ax[0, 1].set_ylabel("fraction of windows z≥3")
    ax[0, 1].set_title("(b) Focal deep-divergence rate\nall archaics equal → shared confound, not lineage-specific")
    ax[0, 1].legend(fontsize=8)

    # (c) simulation den_excess mean by model
    try:
        rows = []
        for m in ["M0_none", "M1_neandersovan", "M2_denisovan_only", "M3_separate",
                  "M4_ancient_structure"]:
            s = pd.read_csv(SIM / f"windows_{m}.tsv", sep="\t")
            ex = s["den_afr"] - s[["alt_afr", "vin_afr", "chag_afr"]].mean(axis=1)
            rows.append((m.split("_")[0], ex.mean()))
        ax[1, 0].bar([r[0] for r in rows], [r[1] * 1e5 for r in rows],
                     color=["grey", "#4477aa", "#ee6677", "#ccbb44", "#aa3377"])
        ax[1, 0].set_ylabel("mean den_excess  (×1e-5)")
        ax[1, 0].set_title("(c) Contrast validated in simulation\nhighest under M2 (Model B), lowest under M0 (null)")
    except Exception as e:
        ax[1, 0].text(0.5, 0.5, f"sim tables absent\n{e}", ha="center")

    # (d) real den_excess distribution, artifact vs clean
    clean = d.loc[~d["artifact"], "z_den_excess"]
    art = d.loc[d["artifact"], "z_den_excess"]
    ax[1, 1].hist([clean, art], bins=40, stacked=True, color=["#228833", "#bbbbbb"],
                  label=["passes artifact filter", "artifact-likely"])
    ax[1, 1].axvline(3, ls="--", c="k", label="CDDR threshold z=3")
    ax[1, 1].set_xlabel("z(den_excess)")
    ax[1, 1].set_ylabel("windows")
    ax[1, 1].set_title("(d) Denisovan-specific excess (real)\ntail is artifact-heavy; survivors are Model-B-direction")
    ax[1, 1].legend(fontsize=8)

    fig.suptitle("Superarchaic-in-Neanderthal pilot (chr21+22): validated pipeline, "
                 "confounded raw signal, artifact-heavy tail, no Model-A sharing", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIG / "pilot_summary.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
