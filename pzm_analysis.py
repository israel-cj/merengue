#!/usr/bin/env python3
"""
pzm_analysis.py - Batch extraction of S-E / P-E loop parameters from aixACCT
TF Analyzer "PiezoResult" (.dat) files (PZM module).

For every measurement table (= one sine period at a given amplitude) it extracts:

  Strain loop (S-E, main figure)
    Smax+  strain at the positive field maximum  (V = V+max)
    Smax-  lowest strain on the negative-field side (bottom of the butterfly
           dip); the strain exactly at V-max is kept as S_at_Vmax-
    Sr     remanent strain: strain where V crosses 0 on the way down
           (after V+max, before V-max)
    Spp    peak-to-peak strain (max - min over the whole loop)

  Polarization loop (P-E, inset)
    Ps+ / Ps-   polarization at V+max / V-max
    Pr+ / Pr-   polarization at V = 0 (descending / ascending branch)
    2Pr         Pr+ - Pr-
    Ec+ / Ec-   field where P crosses 0 (ascending / descending branch)

Thickness and electrode area are read from each table's own header, so
nothing is hard-coded. E [kV/cm] = V / thickness ; S [%] = D / thickness * 100.

Usage
-----
    python pzm_analysis.py  INPUT_FOLDER  OUTPUT_FOLDER  [options]

    options:
      --strain-channel D1     displacement column to use (D1, D2 or D3)
      --pol-channel P1        polarization column to use (P1, P2 or P3)
      --thickness-um 900      override the thickness from the file (um)
      --avg-points 0          average +-N points around V+max / V-max
      --no-plots              only write the CSV / Excel files
      --keep-broken           keep measurements with Measurement Status != 0

Measurements with 'Measurement Status' != 0 are treated as broken: they are
removed from the results and plots and written to broken_measurements*.csv
(use --keep-broken to keep them).

Outputs (in OUTPUT_FOLDER)
-------
    S_vs_X_all_files.csv / .png      X (applied field), Sr, Smax+, Smax- per measurement
                                     (+ S_sum = |Smax+|+|Smax-| and Sr_ratio = Sr/S_sum)
    Sr_ratio_all_files.png           Sr_ratio vs X, and Sr vs S_sum with linear fits
    Sr_ratio_linear_fits.csv         slope, intercept, R2 of Sr vs S_sum per file
    broken_measurements_all_files.csv  measurements with Measurement Status != 0
    summary_all_files.csv / .xlsx    one row per table, all files together
    <file>/<file>_results.csv        same, for one file
    <file>/<file>_loops.csv          all loop points (long format, for replotting)
    <file>/<file>_SE.png, _PE.png    all loops coloured by applied field
    <file>/<file>_params_vs_E.png    Smax+, Smax-, Sr, 2Pr, Ec vs applied field
    <file>/<file>_annotated_max.png  largest loop with Smax+/Smax-/Sr marked (check)

Requires: numpy, pandas, matplotlib (openpyxl optional, for the .xlsx).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# 1. Parsing
# --------------------------------------------------------------------------

TABLE_RE = re.compile(r"^Table\s+(\d+)\s*$")


def _to_float(text: str):
    try:
        return float(text.split()[0])
    except (ValueError, IndexError):
        return text.strip()


def parse_pzm_file(path: Path):
    """Return (summary_df, tables) where tables is a list of dicts:
       {'table': int, 'meta': {...}, 'data': DataFrame}."""
    lines = Path(path).read_text(encoding="latin-1").splitlines()
    n = len(lines)

    # --- summary table (first 'Table 1' followed by a 'Table No' header) ---
    summary = None
    i = 0
    while i < n:
        if TABLE_RE.match(lines[i].strip()) and i + 1 < n and lines[i + 1].startswith("Table No"):
            header = lines[i + 1].rstrip("\t\n ").split("\t")
            rows = []
            j = i + 2
            while j < n and lines[j].strip():
                rows.append(lines[j].rstrip("\t ").split("\t"))
                j += 1
            summary = pd.DataFrame(rows, columns=header).apply(pd.to_numeric, errors="coerce")
            summary = summary.rename(columns={"Table No [#]": "table"})
            i = j
            break
        i += 1

    # --- start of the per-measurement section ('Piezo' line) ---
    while i < n and lines[i].strip() != "Piezo":
        i += 1

    # --- per-measurement tables ---
    tables = []
    while i < n:
        m = TABLE_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue
        tnum = int(m.group(1))
        i += 1
        meta = {}
        # key: value lines until the tab-separated data header
        while i < n and "\t" not in lines[i]:
            line = lines[i].strip()
            if TABLE_RE.match(line):
                break
            if ":" in line:
                k, v = line.split(":", 1)
                meta[k.strip()] = _to_float(v.strip()) if v.strip() else ""
            i += 1
        if i >= n or "\t" not in lines[i]:
            continue  # metadata-only block (e.g. software header)
        header = lines[i].rstrip("\t ").split("\t")
        i += 1
        rows = []
        while i < n and lines[i].strip() and not TABLE_RE.match(lines[i].strip()):
            rows.append(lines[i].rstrip("\t ").split("\t"))
            i += 1
        data = pd.DataFrame(rows, columns=header).apply(pd.to_numeric, errors="coerce")
        tables.append({"table": tnum, "meta": meta, "data": data})

    return summary, tables


# --------------------------------------------------------------------------
# 2. Loop analysis
# --------------------------------------------------------------------------

def _zero_cross(x, y, i0, i1, direction):
    """Linear interpolation of y where x crosses 0 within indices [i0, i1]
    (cyclic: indices wrap around). direction = -1 (falling) or +1 (rising).
    Returns (y_at_crossing, index_float) or (nan, nan)."""
    N = len(x)
    idx = [(i0 + k) % N for k in range(((i1 - i0) % N) + 1)]
    for a, b in zip(idx[:-1], idx[1:]):
        xa, xb = x[a], x[b]
        if direction < 0 and xa > 0 >= xb or direction > 0 and xa < 0 <= xb:
            f = xa / (xa - xb) if xa != xb else 0.0
            return y[a] + f * (y[b] - y[a]), a + f
    return np.nan, np.nan


def _value_at_extreme(v, y, idx, avg):
    lo, hi = max(idx - avg, 0), min(idx + avg + 1, len(v))
    return float(np.mean(y[lo:hi])), float(np.mean(v[lo:hi]))


def analyse_table(t, strain_ch="D1", pol_ch="P1", thickness_um=None, avg=0):
    meta, d = t["meta"], t["data"]
    th_nm = thickness_um * 1e3 if thickness_um else float(meta.get("Thickness [nm]", np.nan))
    th_cm = th_nm * 1e-7

    V = d["V+ [V]"].to_numpy(float)
    D = d[f"{strain_ch} [nm]"].to_numpy(float)
    P = d[f"{pol_ch} [uC/cm2]"].to_numpy(float)
    E = V / th_cm / 1e3                 # kV/cm
    S = D / th_nm * 100.0               # %

    imax, imin = int(np.argmax(V)), int(np.argmin(V))
    Smax_p, Emax_p = _value_at_extreme(E, S, imax, avg)
    S_at_Vmin, Emax_m = _value_at_extreme(E, S, imin, avg)
    # Smax-: lowest strain on the negative-field side (bottom of the butterfly
    # dip). At low field this coincides with V-max; after switching it does not.
    neg = np.where(E < 0)[0]
    ineg = int(neg[np.nanargmin(S[neg])]) if len(neg) else imin
    Smax_m, E_Smax_m = _value_at_extreme(E, S, ineg, avg)
    Ps_p, _ = _value_at_extreme(E, P, imax, avg)
    Ps_m, _ = _value_at_extreme(E, P, imin, avg)

    # descending branch: V+max -> V-max ; ascending: V-max -> V+max (wraps)
    Sr, _ = _zero_cross(V, S, imax, imin, -1)
    Pr_p, _ = _zero_cross(V, P, imax, imin, -1)
    Pr_m, _ = _zero_cross(V, P, imin, imax, +1)
    S_end, _ = _zero_cross(V, S, imin, imax, +1)     # strain back at V=0
    Ec_m, _ = _zero_cross(P, E, imax, imin, -1)      # P falls through 0
    Ec_p, _ = _zero_cross(P, E, imin, imax, +1)      # P rises through 0

    return {
        "table": t["table"],
        "sample_name": meta.get("SampleName", ""),
        "timestamp": meta.get("Timestamp", ""),
        "waveform": meta.get("Waveform", ""),
        "freq_Hz": meta.get("Hysteresis Frequency [Hz]", np.nan),
        "amplitude_V": meta.get("Hysteresis Amplitude [V]", np.nan),
        "thickness_um": th_nm / 1e3,
        "area_mm2": meta.get("Area [mm2]", np.nan),
        "Eapp_kV_cm": meta.get("Hysteresis Amplitude [V]", np.nan) / th_cm / 1e3,
        "Emax+_kV_cm": Emax_p,
        "Emax-_kV_cm": Emax_m,
        "Smax+_%": Smax_p,
        "Smax-_%": Smax_m,
        "E_at_Smax-_kV_cm": E_Smax_m,
        "S_at_Vmax-_%": S_at_Vmin,
        "Sr_%": Sr,
        "S_end_%": S_end,
        "Spp_%": float(np.nanmax(S) - np.nanmin(S)),
        # true minimum of the loop (butterfly dip) - differs from Smax- once
        # the sample switches at high field
        "Smin_%": float(np.nanmin(S)),
        "E_at_Smin_kV_cm": float(E[np.nanargmin(S)]),
        "Ps+_uC_cm2": Ps_p,
        "Ps-_uC_cm2": Ps_m,
        "Pr+_uC_cm2": Pr_p,
        "Pr-_uC_cm2": Pr_m,
        "2Pr_uC_cm2": Pr_p - Pr_m,
        "Ec+_kV_cm": Ec_p,
        "Ec-_kV_cm": Ec_m,
        "d33*_pm_V": (Smax_p / 100) / (Emax_p * 1e5) * 1e12 if Emax_p else np.nan,
        # instrument values, for cross-checking
        "instr_Dvmax+_nm": meta.get("Dvmax+ [nm/V]", np.nan),
        "instr_Dvmax-_nm": meta.get("Dvmax- [nm/V]", np.nan),
        "instr_Pr+_uC_cm2": meta.get("Pr+ [uC/cm2]", np.nan),
        "instr_Vc+_V": meta.get("Vc+ [V]", np.nan),
        # 0 = OK; any other value (1, 2, ...) = broken measurement
        "measurement_status": meta.get("Measurement Status", np.nan),
    }, pd.DataFrame({"table": t["table"], "t_s": d["Time [s]"], "V": V,
                     "E_kV_cm": E, "D_nm": D, "S_%": S, "P_uC_cm2": P})


def label_sequence(res: pd.DataFrame) -> pd.DataFrame:
    """Mark each table as ramp 'up' or 'down' and flag repeated amplitudes."""
    amp = res["amplitude_V"].to_numpy()
    direction = ["up"]
    for k in range(1, len(amp)):
        if amp[k] > amp[k - 1]:
            direction.append("up")
        elif amp[k] < amp[k - 1]:
            direction.append("down")
        else:
            direction.append("repeat")
    res.insert(1, "step", range(1, len(res) + 1))
    res.insert(2, "ramp", direction)
    # branch: 'ascending' up to (and including) the first loop at the highest
    # amplitude, 'descending' after it. Used for plotting against X.
    ipeak = int(np.nanargmax(amp))
    res.insert(3, "branch", ["ascending" if k <= ipeak else "descending"
                             for k in range(len(amp))])
    # when the same amplitude is measured several times in a row, only the last
    # one (the most stabilised loop) gets is_last_at_amplitude = True
    nxt = np.r_[amp[1:], np.nan]
    res["is_last_at_amplitude"] = amp != nxt
    return res


def is_ok(res: pd.DataFrame) -> pd.Series:
    """True for good measurements: Measurement Status == 0 (or not reported)."""
    st = pd.to_numeric(res["measurement_status"], errors="coerce")
    return st.isna() | (st == 0)


# columns of the compact per-measurement table (X, Sr, Smax+, Smax-)
S_VS_X_COLS = ["file", "step", "table", "branch", "ramp", "amplitude_V",
               "X_kV_cm", "Sr_%", "Smax+_%", "Smax-_%", "E_at_Smax-_kV_cm",
               "S_sum_%", "Sr_ratio", "measurement_status"]


# --------------------------------------------------------------------------
# 3. Plots
# --------------------------------------------------------------------------

def make_plots(stem, res, loops, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    norm = colors.Normalize(res["Eapp_kV_cm"].min(), res["Eapp_kV_cm"].max())
    cmap = cm.viridis
    eapp = dict(zip(res["table"], res["Eapp_kV_cm"]))

    for ycol, ylabel, tag in [("S_%", "Strain (%)", "SE"),
                              ("P_uC_cm2", r"Polarization ($\mu$C/cm$^2$)", "PE")]:
        fig, ax = plt.subplots(figsize=(6, 5))
        for tnum, g in loops.groupby("table"):
            x, y = np.r_[g["E_kV_cm"], g["E_kV_cm"].iloc[0]], np.r_[g[ycol], g[ycol].iloc[0]]
            ax.plot(x, y, lw=1, color=cmap(norm(eapp[tnum])))
        ax.axhline(0, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6)
        ax.set_xlabel("E (kV/cm)"); ax.set_ylabel(ylabel); ax.set_title(stem, fontsize=9)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="E$_{app}$ (kV/cm)")
        fig.tight_layout(); fig.savefig(outdir / f"{stem}_{tag}.png", dpi=200); plt.close(fig)

    # parameters vs applied field (up / down ramps distinguished)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    marker = {"up": "o-", "down": "s--", "repeat": "x"}
    for ramp, g in res.groupby("ramp"):
        m = marker.get(ramp, "o")
        axs[0].plot(g["Eapp_kV_cm"], g["Smax+_%"], m, label=f"Smax+ ({ramp})", color="C0")
        axs[0].plot(g["Eapp_kV_cm"], g["Sr_%"], m, label=f"Sr ({ramp})", color="C1")
        axs[0].plot(g["Eapp_kV_cm"], g["Smax-_%"], m, label=f"Smax- ({ramp})", color="C2")
        axs[1].plot(g["Eapp_kV_cm"], g["2Pr_uC_cm2"], m, label=f"2Pr ({ramp})", color="C3")
        axs[2].plot(g["Eapp_kV_cm"], g["Ec+_kV_cm"], m, label=f"Ec+ ({ramp})", color="C4")
        axs[2].plot(g["Eapp_kV_cm"], g["Ec-_kV_cm"], m, label=f"Ec- ({ramp})", color="C5")
    for ax, yl in zip(axs, ["Strain (%)", r"2P$_r$ ($\mu$C/cm$^2$)", "E$_c$ (kV/cm)"]):
        ax.set_xlabel("E$_{app}$ (kV/cm)"); ax.set_ylabel(yl); ax.legend(fontsize=7)
        ax.axhline(0, color="grey", lw=0.5)
    fig.suptitle(stem, fontsize=9); fig.tight_layout()
    fig.savefig(outdir / f"{stem}_params_vs_E.png", dpi=200); plt.close(fig)

    # annotated largest loop: visual check that the points are picked correctly
    r = res.loc[res["amplitude_V"].idxmax()]
    g = loops[loops["table"] == r["table"]]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(np.r_[g["E_kV_cm"], g["E_kV_cm"].iloc[0]], np.r_[g["S_%"], g["S_%"].iloc[0]], color="navy")
    ax.axhline(0, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6)
    pts = [(r["Emax+_kV_cm"], r["Smax+_%"], "S$_{max}^+$"),
           (r["E_at_Smax-_kV_cm"], r["Smax-_%"], "S$_{max}^-$"),
           (0.0, r["Sr_%"], "S$_r$")]
    for x, y, lab in pts:
        ax.plot(x, y, "o", color="crimson", zorder=5)
        ax.annotate(f"{lab} = {y:.4f} %", (x, y), textcoords="offset points",
                    xytext=(8, 8), fontsize=9, color="crimson")
    ax.set_xlabel("E (kV/cm)"); ax.set_ylabel("Strain (%)")
    ax.set_title(f"{stem} - table {int(r['table'])}, E$_{{app}}$ = {r['Eapp_kV_cm']:.1f} kV/cm", fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / f"{stem}_annotated_max.png", dpi=200); plt.close(fig)


SERIES = [  # (column, label, colour, marker) - fixed order, colour-blind safe
    ("Smax+_%", r"S$_{max}^{+}$", "#2a78d6", "o"),
    ("Sr_%",    r"S$_{r}$",       "#eb6834", "s"),
    ("Smax-_%", r"S$_{max}^{-}$", "#1baf7a", "^"),
]


def plot_s_vs_x(df, title, outpath, n_broken=0):
    """Sr, Smax+, Smax- vs X (applied field) for one file.
    When an amplitude is measured several times in a row only the last loop is
    plotted (is_last_at_amplitude); all of them stay in the CSV files.
    Ascending branch: filled markers + solid line; descending: open + dashed."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # repeated amplitudes: plot only the last (stabilised) loop; CSVs keep all
    df = df[df["is_last_at_amplitude"]]
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for col, lab, colr, mk in SERIES:
        for br, ls, face in [("ascending", "-", colr), ("descending", "--", "white")]:
            g = df[df["branch"] == br].sort_values("step")
            if g.empty:
                continue
            ax.plot(g["X_kV_cm"], g[col], ls=ls, lw=1.5, color=colr, marker=mk,
                    ms=6, mfc=face, mec=colr, mew=1.5,
                    label=f"{lab} ({br})" if df["branch"].nunique() > 1 else lab)
        # direct label at the end of the ascending branch
        g = df[df["branch"] == "ascending"].sort_values("step")
        if not g.empty:
            ax.annotate(lab, (g["X_kV_cm"].iloc[-1], g[col].iloc[-1]),
                        textcoords="offset points", xytext=(6, 0), va="center",
                        fontsize=10, color="#333333")
    ax.axhline(0, color="#999999", lw=0.8)
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("X = applied electric field E$_{app}$ (kV/cm)")
    ax.set_ylabel("Strain (%)")
    note = f"  ({n_broken} broken measurement(s) excluded)" if n_broken else ""
    ax.set_title(title + note, fontsize=9)
    ax.legend(fontsize=8, frameon=False, ncol=2 if df["branch"].nunique() > 1 else 1)
    fig.tight_layout(); fig.savefig(outpath, dpi=200); plt.close(fig)


def plot_s_vs_x_all_files(df, outpath):
    """One panel per quantity (Smax+, Sr, Smax-), one line per file."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = df[df["is_last_at_amplitude"]]   # see plot_s_vs_x
    files = list(dict.fromkeys(df["file"]))
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
    markers = ["o", "s", "^", "D", "v", "P"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), sharex=True)
    for ax, (col, lab, _, _) in zip(axs, SERIES):
        for k, f in enumerate(files):
            c, m = colors[k % len(colors)], markers[k % len(markers)]
            g = df[df["file"] == f]
            for br, ls, face in [("ascending", "-", c), ("descending", "--", "white")]:
                h = g[g["branch"] == br].sort_values("step")
                if not h.empty:
                    ax.plot(h["X_kV_cm"], h[col], ls=ls, lw=1.3, color=c, marker=m,
                            ms=5, mfc=face, mec=c,
                            label=Path(f).stem if br == "ascending" else None)
        ax.axhline(0, color="#999999", lw=0.8)
        ax.grid(True, color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_title(lab, fontsize=12)
        ax.set_xlabel("X = E$_{app}$ (kV/cm)")
    axs[0].set_ylabel("Strain (%)")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center",
               ncol=min(len(files), 3), fontsize=8, frameon=False)
    fig.text(0.99, 0.01, "filled/solid = ascending field, open/dashed = descending",
             ha="right", fontsize=7, color="#555555")
    fig.tight_layout(rect=(0, 0.1 + 0.04 * ((len(files) - 1) // 3), 1, 1))
    fig.savefig(outpath, dpi=200); plt.close(fig)


FILE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
FILE_MARKERS = ["o", "s", "^", "D", "v", "P"]


def linear_fits(df):
    """Least-squares line Sr = slope * (|Smax+| + |Smax-|) + intercept, per file.
    Uses the same points as the plots (last loop at each amplitude)."""
    rows = []
    for f, g in df[df["is_last_at_amplitude"]].groupby("file", sort=False):
        x, y = g["S_sum_%"].to_numpy(float), g["Sr_%"].to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y)
        x, y = x[m], y[m]
        if len(x) < 3:
            continue
        slope, icpt = np.polyfit(x, y, 1)
        r2 = 1 - np.sum((y - (slope * x + icpt)) ** 2) / np.sum((y - y.mean()) ** 2)
        rows.append({"file": f, "n_points": len(x), "slope": slope,
                     "intercept_%": icpt, "R2": r2,
                     "Sr_ratio_mean": float(np.mean(y / x)),
                     "Sr_ratio_std": float(np.std(y / x, ddof=1))})
    return pd.DataFrame(rows)


def plot_sr_ratio(df, fits, title, outpath):
    """(a) Sr / (|Smax+| + |Smax-|) vs X.
       (b) Sr vs (|Smax+| + |Smax-|) with a linear fit per file.
    One colour per file; filled/solid = ascending, open/dashed = descending."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = df[df["is_last_at_amplitude"]]
    files = list(dict.fromkeys(df["file"]))
    single = len(files) == 1
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.9))
    for k, f in enumerate(files):
        c, mk = FILE_COLORS[k % len(FILE_COLORS)], FILE_MARKERS[k % len(FILE_MARKERS)]
        g = df[df["file"] == f]
        for br, ls, face in [("ascending", "-", c), ("descending", "--", "white")]:
            h = g[g["branch"] == br].sort_values("step")
            if h.empty:
                continue
            lab = (br if single else Path(f).stem) if (single or br == "ascending") else None
            a1.plot(h["X_kV_cm"], h["Sr_ratio"], ls=ls, lw=1.3, color=c, marker=mk,
                    ms=6, mfc=face, mec=c, mew=1.3, label=lab)
            a2.plot(h["S_sum_%"], h["Sr_%"], ls="none", color=c, marker=mk,
                    ms=6, mfc=face, mec=c, mew=1.3)
        fr = fits[fits["file"] == f]
        if not fr.empty:
            r = fr.iloc[0]
            xx = np.linspace(0, g["S_sum_%"].max() * 1.05, 50)
            a2.plot(xx, r["slope"] * xx + r["intercept_%"], color=c, lw=1.2,
                    label=(f"{'' if single else Path(f).stem + ': '}"
                           f"slope {r['slope']:.3f}, R$^2$ = {r['R2']:.3f}"))
    a1.set_xlabel("X = applied electric field E$_{app}$ (kV/cm)")
    a1.set_ylabel(r"S$_r$ / (|S$_{max}^{+}$| + |S$_{max}^{-}$|)")
    a1.set_title(r"(a) Normalised remanent strain vs X", fontsize=10)
    a2.set_xlabel(r"|S$_{max}^{+}$| + |S$_{max}^{-}$|  (%)")
    a2.set_ylabel(r"S$_r$ (%)")
    a2.set_title("(b) S$_r$ vs total strain, linear fit", fontsize=10)
    a2.set_xlim(left=0); a2.set_ylim(bottom=min(0, df["Sr_%"].min() * 1.1))
    for ax in (a1, a2):
        ax.grid(True, color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle(title + "   (filled/solid = ascending field, open/dashed = descending)",
                 fontsize=9)
    fig.tight_layout(); fig.savefig(outpath, dpi=200); plt.close(fig)


# --------------------------------------------------------------------------
# 4. Driver
# --------------------------------------------------------------------------

def process_file(path, outroot, args):
    stem = path.stem
    summary, tables = parse_pzm_file(path)
    if not tables:
        print(f"  ! {path.name}: no measurement tables found, skipped")
        return None
    rows, loops = [], []
    for t in tables:
        try:
            r, lp = analyse_table(t, args.strain_channel, args.pol_channel,
                                  args.thickness_um, args.avg_points)
            rows.append(r); loops.append(lp)
        except Exception as e:  # keep going if one table is broken
            print(f"  ! {path.name} table {t['table']}: {e}")
    res = label_sequence(pd.DataFrame(rows))
    res.insert(0, "file", path.name)
    res["X_kV_cm"] = res["Eapp_kV_cm"]          # X = applied electric field
    # Sr normalised by the total strain swing |Smax+| + |Smax-|
    res["S_sum_%"] = res["Smax+_%"].abs() + res["Smax-_%"].abs()
    res["Sr_ratio"] = res["Sr_%"] / res["S_sum_%"]
    loops = pd.concat(loops, ignore_index=True)

    # --- split good / broken (Measurement Status != 0) ---
    ok = is_ok(res)
    broken = res[~ok].copy()
    good = res if args.keep_broken else res[ok].copy()
    good_loops = loops[loops["table"].isin(good["table"])]
    broken_loops = loops[loops["table"].isin(broken["table"])]
    if not args.keep_broken and len(good):
        # recompute "last loop at this amplitude" without the broken ones
        a = good["amplitude_V"].to_numpy()
        good["is_last_at_amplitude"] = a != np.r_[a[1:], np.nan]

    outdir = outroot / stem
    outdir.mkdir(parents=True, exist_ok=True)
    good.to_csv(outdir / f"{stem}_results.csv", index=False)
    good[S_VS_X_COLS].to_csv(outdir / f"{stem}_S_vs_X.csv", index=False)
    good_loops.to_csv(outdir / f"{stem}_loops.csv", index=False)
    broken.to_csv(outdir / f"{stem}_broken_measurements.csv", index=False)
    if len(broken_loops):
        broken_loops.to_csv(outdir / f"{stem}_broken_loops.csv", index=False)
    if summary is not None:
        summary.to_csv(outdir / f"{stem}_instrument_summary.csv", index=False)
    if not args.no_plots and len(good):
        make_plots(stem, good, good_loops, outdir)
        plot_s_vs_x(good, stem, outdir / f"{stem}_S_vs_X.png",
                    n_broken=0 if args.keep_broken else len(broken))
        plot_sr_ratio(good, linear_fits(good), stem, outdir / f"{stem}_Sr_ratio.png")
    msg = f", {len(broken)} broken (status != 0)" if len(broken) else ""
    print(f"  {path.name}: {len(res)} loops{msg}, thickness {res['thickness_um'].iloc[0]:.0f} um, "
          f"Eapp {res['Eapp_kV_cm'].min():.1f}-{res['Eapp_kV_cm'].max():.1f} kV/cm")
    return good, broken


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder", type=Path)
    ap.add_argument("output_folder", type=Path)
    ap.add_argument("--strain-channel", default="D1", choices=["D1", "D2", "D3"])
    ap.add_argument("--pol-channel", default="P1", choices=["P1", "P2", "P3"])
    ap.add_argument("--thickness-um", type=float, default=None)
    ap.add_argument("--avg-points", type=int, default=0)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--keep-broken", action="store_true",
                    help="keep measurements with Measurement Status != 0 in the results "
                         "(they are still listed in the broken_measurements CSV)")
    args = ap.parse_args(argv)

    files = sorted(p for p in args.input_folder.iterdir() if p.suffix.lower() == ".dat")
    if not files:
        sys.exit(f"No .dat files in {args.input_folder}")
    args.output_folder.mkdir(parents=True, exist_ok=True)
    print(f"Processing {len(files)} file(s) -> {args.output_folder}")

    out = [r for f in files if (r := process_file(f, args.output_folder, args)) is not None]
    allres = pd.concat([g for g, _ in out], ignore_index=True)
    allbroken = pd.concat([b for _, b in out], ignore_index=True)
    s_vs_x = allres[S_VS_X_COLS]
    o = args.output_folder
    allres.to_csv(o / "summary_all_files.csv", index=False)
    s_vs_x.to_csv(o / "S_vs_X_all_files.csv", index=False)
    allbroken.to_csv(o / "broken_measurements_all_files.csv", index=False)
    fits = linear_fits(allres)
    fits.to_csv(o / "Sr_ratio_linear_fits.csv", index=False)
    if not args.no_plots and len(allres):
        plot_s_vs_x_all_files(allres, o / "S_vs_X_all_files.png")
        plot_sr_ratio(allres, fits, "All files", o / "Sr_ratio_all_files.png")
    try:
        with pd.ExcelWriter(o / "summary_all_files.xlsx") as xw:
            s_vs_x.to_excel(xw, sheet_name="S_vs_X", index=False)
            allres.to_excel(xw, sheet_name="all_parameters", index=False)
            fits.to_excel(xw, sheet_name="Sr_ratio_linear_fits", index=False)
            allbroken.to_excel(xw, sheet_name="broken_measurements", index=False)
            for fname, g in allres.groupby("file"):
                g[S_VS_X_COLS].to_excel(xw, sheet_name=Path(fname).stem[-31:], index=False)
    except ImportError:
        print("  (openpyxl not installed: Excel summary skipped, CSV written)")
    print(f"Done. {len(allres)} measurements kept, {len(allbroken)} broken "
          f"measurement(s) {'flagged' if args.keep_broken else 'excluded'} "
          f"-> broken_measurements_all_files.csv")


if __name__ == "__main__":
    main()
