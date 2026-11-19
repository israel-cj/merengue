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
    python pzm_analysis.py  INPUT_FOLDER_OR_FILE  OUTPUT_FOLDER  [options]

    options:
      --strain-channel D1     displacement column to use (D1, D2 or D3)
      --pol-channel P1        polarization column to use (P1, P2 or P3)
      --thickness-um 900      override the thickness from the file (um)
      --avg-points 0          average +-N points around V+max / V-max
      --no-plots              only write the CSV / Excel files
      --keep-broken           keep measurements with Measurement Status != 0
      --slope-window -10 10   field window (kV/cm) for the slope of D [nm] vs V+ [V]
      --slope-eapp all        which loops go in the slope table, by E_app:
                              20 | 10:20 | >=15 | "10:15, 20" | all
      --interactive           show the E_app values found and ask for the
                              slope window and E_app condition

Measurements with 'Measurement Status' != 0 are treated as broken: they are
removed from the results and plots and written to broken_measurements*.csv
(use --keep-broken to keep them).

Outputs (in OUTPUT_FOLDER)
-------
    S_vs_X_all_files.csv / .png      X (applied field), Sr, Smax+, Smax- per measurement
                                     (+ S_sum = |Smax+|+|Smax-| and Sr_ratio = Sr/S_sum)
    Sr_ratio_all_files.png           Sr_ratio vs X, and Sr vs S_sum with linear fits
    Sr_ratio_linear_fits.csv         slope, intercept, R2 of Sr vs S_sum per file
    slopes_all_files.csv / .png      slope dD/dV (nm/V) in the field window, per loop meeting
                                     the E_app condition (<file>/<file>_slopes.csv,
                                     <file>/<file>_slope_check.png per file)
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


def _linfit(x, y):
    """Least-squares line y = slope * x + icpt. Returns (slope, icpt, R2, n)."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3 or np.ptp(x) == 0:
        return np.nan, np.nan, np.nan, len(x)
    slope, icpt = np.polyfit(x, y, 1)
    ss = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum((y - (slope * x + icpt)) ** 2) / ss if ss > 0 else np.nan
    return float(slope), float(icpt), float(r2), len(x)


def window_slopes(E, x, y, imax, imin, e_lo, e_hi):
    """Slope dy/dx of a linear fit to the points with e_lo <= E <= e_hi (the
    window is given in field, kV/cm), done separately on each half-cycle:
        'down' = decreasing field, V+max -> V-max
        'up'   = increasing field, V-max -> V+max (wraps through t = 0)
    Used with x = V+ [V] and y = D [nm]  ->  slope in nm/V.
    Returns {branch: (slope, icpt, R2, n, indices)}."""
    N = len(E)
    out = {}
    for name, (a, b) in (("down", (imax, imin)), ("up", (imin, imax))):
        idx = np.array([(a + k) % N for k in range(((b - a) % N) + 1)])
        sel = idx[(E[idx] >= e_lo) & (E[idx] <= e_hi)]
        out[name] = _linfit(x[sel], y[sel]) + (sel,)
    return out


# the loop must reach the window edges to within this fraction of the window
# width (E at V+max is a sampled value, e.g. 9.98 kV/cm for E_app = 10 kV/cm)
WINDOW_TOL = 0.01


def analyse_table(t, strain_ch="D1", pol_ch="P1", thickness_um=None, avg=0,
                  slope_window=(-10.0, 10.0)):
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

    # slope dD/dV (raw V+ [V] vs D [nm], nm/V) for the points inside the field
    # window [e_lo, e_hi] kV/cm  (= [e_lo, e_hi] * thickness in volts)
    e_lo, e_hi = slope_window
    ws = window_slopes(E, V, D, imax, imin, e_lo, e_hi)
    tol = WINDOW_TOL * (e_hi - e_lo)
    covered = bool(np.nanmax(E) >= e_hi - tol and np.nanmin(E) <= e_lo + tol)
    sl_up, sl_dn = ws["up"][0], ws["down"][0]
    sl_mean = np.nanmean([sl_up, sl_dn]) if np.isfinite([sl_up, sl_dn]).any() else np.nan

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
        # slope of D [nm] vs V+ [V] for the points inside the field window
        "slope_Elo_kV_cm": e_lo,
        "slope_Ehi_kV_cm": e_hi,
        "slope_Vlo_V": e_lo * 1e3 * th_cm,
        "slope_Vhi_V": e_hi * 1e3 * th_cm,
        "slope_window_covered": covered,
        "slope_up_nm_per_V": sl_up,
        "slope_down_nm_per_V": sl_dn,
        "slope_mean_nm_per_V": sl_mean,
        "slope_up_R2": ws["up"][2],
        "slope_down_R2": ws["down"][2],
        "slope_up_n": ws["up"][3],
        "slope_down_n": ws["down"][3],
        "slope_up_icpt_nm": ws["up"][1],
        "slope_down_icpt_nm": ws["down"][1],
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

# columns of the slope table (dD/dV, V+ [V] vs D [nm], inside the field window)
SLOPE_COLS = ["file", "step", "table", "branch", "ramp", "amplitude_V", "X_kV_cm",
              "thickness_um", "Emax+_kV_cm", "Emax-_kV_cm",
              "slope_Elo_kV_cm", "slope_Ehi_kV_cm", "slope_Vlo_V", "slope_Vhi_V",
              "slope_window_covered", "slope_up_nm_per_V", "slope_down_nm_per_V",
              "slope_mean_nm_per_V", "slope_up_R2", "slope_down_R2",
              "slope_up_n", "slope_down_n", "is_last_at_amplitude", "measurement_status"]

EAPP_TOL = 0.05   # kV/cm, tolerance when matching E_app values (8.0000001 == 8)


def parse_eapp_condition(spec):
    """Turn a text condition on the applied field X = E_app (kV/cm) into a
    function array -> bool mask. Accepted forms (combine with commas):
        all | (empty)     every loop
        20                E_app == 20 (within EAPP_TOL)
        10:20  or 10-20   10 <= E_app <= 20
        >=15, >15, <=12, <12
    e.g. "10:15, 20"  ->  10 <= E_app <= 15  or  E_app == 20"""
    spec = (spec or "").strip().lower()
    if spec in ("", "all", "*"):
        return lambda x: np.ones(len(x), bool)
    tests = []
    for part in (p.strip() for p in spec.split(",")):
        if not part:
            continue
        m = re.fullmatch(r"(>=|<=|>|<)\s*([\d.]+)", part)
        r = re.fullmatch(r"([\d.]+)\s*[:-]\s*([\d.]+)", part)
        try:
            if m:
                op, v = m.group(1), float(m.group(2))
                tests.append({">=": lambda x, v=v: x >= v - EAPP_TOL,
                              "<=": lambda x, v=v: x <= v + EAPP_TOL,
                              ">": lambda x, v=v: x > v + EAPP_TOL,
                              "<": lambda x, v=v: x < v - EAPP_TOL}[op])
            elif r:
                lo, hi = sorted(map(float, r.groups()))
                tests.append(lambda x, lo=lo, hi=hi: (x >= lo - EAPP_TOL) & (x <= hi + EAPP_TOL))
            else:
                v = float(part)
                tests.append(lambda x, v=v: np.abs(x - v) <= EAPP_TOL)
        except ValueError:
            raise ValueError(f"cannot read E_app condition '{part}' "
                             "(use e.g. 20, 10:20, >=15, or a comma list)") from None
    return lambda x: np.logical_or.reduce([t(np.asarray(x, float)) for t in tests])


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


FILE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
               "#7a4fb5", "#4d4d4d"]
FILE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]


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


def plot_slope_check(stem, row, loop, outpath=None, d_label="D1"):
    """One raw D-V loop with the window shaded and the two fitted lines drawn
    on top: visual check of the slope calculation. Returns the figure when
    outpath is None (used by the Gradio app), otherwise saves it."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    v_lo, v_hi = row["slope_Vlo_V"], row["slope_Vhi_V"]
    V, D = loop["V"].to_numpy(), loop["D_nm"].to_numpy()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.axvspan(v_lo, v_hi, color="#2a78d6", alpha=0.08, lw=0)
    ax.plot(np.r_[V, V[0]], np.r_[D, D[0]], color="#888888", lw=1, label="D-V loop")
    xx = np.array([v_lo, v_hi])
    for br, colr, lab in [("up", "#eb6834", "increasing V"), ("down", "#1baf7a", "decreasing V")]:
        m, b = row[f"slope_{br}_nm_per_V"], row[f"slope_{br}_icpt_nm"]
        if np.isfinite(m):
            ax.plot(xx, m * xx + b, color=colr, lw=2,
                    label=f"{lab}: {m:.4f} nm/V, R$^2$ = {row[f'slope_{br}_R2']:.3f}")
    ax.axhline(0, color="#bbbbbb", lw=0.6); ax.axvline(0, color="#bbbbbb", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("V+ (V)"); ax.set_ylabel(f"{d_label} (nm)")
    ax.set_title(f"{stem} - table {int(row['table'])}, E$_{{app}}$ = {row['X_kV_cm']:.1f} kV/cm\n"
                 f"fit in {row['slope_Elo_kV_cm']:g} ... {row['slope_Ehi_kV_cm']:g} kV/cm "
                 f"= {v_lo:.0f} ... {v_hi:.0f} V (shaded), "
                 f"mean {row['slope_mean_nm_per_V']:.4f} nm/V", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.tight_layout()
    if outpath is None:
        return fig
    fig.savefig(outpath, dpi=200); plt.close(fig)


def plot_slopes_all_files(df, outpath, e_lo, e_hi, condition):
    """Mean slope dD/dV (nm/V) vs X, one colour per file. Points are joined in
    measurement order; filled/solid = amplitude ramping up, open/dashed = down.
    Returns the figure when outpath is None."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = df[df["is_last_at_amplitude"] & df["slope_window_covered"]]
    if df.empty:
        return None
    files = list(dict.fromkeys(df["file"]))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for k, f in enumerate(files):
        c, mk = FILE_COLORS[k % len(FILE_COLORS)], FILE_MARKERS[k % len(FILE_MARKERS)]
        g = df[df["file"] == f].sort_values("step")
        down = (g["ramp"] == "down").to_numpy()
        run = np.r_[0, np.cumsum(down[1:] != down[:-1])]   # consecutive up / down runs
        lab = Path(f).stem                                   # legend entry once per file
        for k_run in np.unique(run):
            idx = np.flatnonzero(run == k_run)
            is_down = down[idx[0]]
            line = np.r_[idx[0] - 1, idx] if idx[0] > 0 else idx   # join to previous run
            ax.plot(g["X_kV_cm"].iloc[line], g["slope_mean_nm_per_V"].iloc[line],
                    ls="--" if is_down else "-", lw=1.3, color=c, zorder=1)
            ax.plot(g["X_kV_cm"].iloc[idx], g["slope_mean_nm_per_V"].iloc[idx], ls="none",
                    marker=mk, ms=6, mfc="white" if is_down else c, mec=c, mew=1.3,
                    zorder=2, label=lab)
            lab = None
    ax.grid(True, color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("X = applied electric field E$_{app}$ (kV/cm)")
    ax.set_ylabel(f"dD/dV for E in [{e_lo:g}, {e_hi:g}] kV/cm  (nm/V)")
    ax.set_title(f"Slope of D vs V+ between {e_lo:g} and {e_hi:g} kV/cm "
                 f"(E$_{{app}}$ condition: {condition or 'all'})\n"
                 "mean of both half-cycles; filled/solid = E$_{app}$ ramping up, "
                 "open/dashed = down", fontsize=9)
    ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    if outpath is None:
        return fig
    fig.savefig(outpath, dpi=200); plt.close(fig)


# --------------------------------------------------------------------------
# 4. Driver
# --------------------------------------------------------------------------

def analyse_file(path, strain_ch="D1", pol_ch="P1", thickness_um=None, avg=0,
                 slope_window=(-10.0, 10.0), keep_broken=False, log=print):
    """Parse and analyse one .dat file, without writing anything.
    Returns (res, good, broken, loops, summary) or None if no tables.
      res    all loops        good   loops kept (status == 0 unless keep_broken)
      broken status != 0      loops  all loop points (long format)"""
    path = Path(path)
    summary, tables = parse_pzm_file(path)
    if not tables:
        return None
    rows, loops = [], []
    for t in tables:
        try:
            r, lp = analyse_table(t, strain_ch, pol_ch, thickness_um, avg,
                                  slope_window=tuple(slope_window))
            rows.append(r); loops.append(lp)
        except Exception as e:  # keep going if one table is broken
            log(f"  ! {path.name} table {t['table']}: {e}")
    if not rows:
        return None
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
    good = res.copy() if keep_broken else res[ok].copy()
    if not keep_broken and len(good):
        # recompute "last loop at this amplitude" without the broken ones
        a = good["amplitude_V"].to_numpy()
        good["is_last_at_amplitude"] = a != np.r_[a[1:], np.nan]
    return res, good, broken, loops, summary


def process_file(path, outroot, args):
    stem = path.stem
    out = analyse_file(path, args.strain_channel, args.pol_channel, args.thickness_um,
                       args.avg_points, tuple(args.slope_window), args.keep_broken)
    if out is None:
        print(f"  ! {path.name}: no measurement tables found, skipped")
        return None
    res, good, broken, loops, summary = out
    good_loops = loops[loops["table"].isin(good["table"])]
    broken_loops = loops[loops["table"].isin(broken["table"])]

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
    # slope table for the loops that meet the E_app condition
    sel = good[args.eapp_test(good["X_kV_cm"].to_numpy())]
    sel[SLOPE_COLS].to_csv(outdir / f"{stem}_slopes.csv", index=False)
    usable = sel[sel["slope_window_covered"]]
    if not args.no_plots and len(usable):
        # check plot: the largest loop that meets the condition
        r = usable.loc[usable["amplitude_V"].idxmax()]
        plot_slope_check(stem, r, good_loops[good_loops["table"] == r["table"]],
                         outdir / f"{stem}_slope_check.png", d_label=args.strain_channel)
    msg = f", {len(broken)} broken (status != 0)" if len(broken) else ""
    print(f"  {path.name}: {len(res)} loops{msg}, thickness {res['thickness_um'].iloc[0]:.0f} um, "
          f"Eapp {res['Eapp_kV_cm'].min():.1f}-{res['Eapp_kV_cm'].max():.1f} kV/cm")
    return good, broken


def _ask(prompt, default):
    """input() with a default; falls back to the default when there is no
    console (e.g. stdin closed)."""
    try:
        ans = input(f"{prompt} [{default}]: ").strip().lstrip("\ufeff")
    except EOFError:
        ans = ""
    return ans or str(default)


def ask_slope_settings(files, args):
    """--interactive: list the E_app values in the files, then ask for the
    slope window and the E_app condition."""
    print("\nApplied fields (E_app, kV/cm) found in the files:")
    for f in files:
        _, tables = parse_pzm_file(f)
        vals = []
        for t in tables:
            th = args.thickness_um * 1e3 if args.thickness_um else t["meta"].get("Thickness [nm]")
            amp = t["meta"].get("Hysteresis Amplitude [V]")
            if isinstance(th, float) and isinstance(amp, float) and th:
                vals.append(round(amp / (th * 1e-7) / 1e3, 1))
        print(f"  {f.stem}: {', '.join(f'{v:g}' for v in sorted(set(vals)))}")
    lo_d, hi_d = args.slope_window
    while True:
        try:
            lo = float(_ask("\nLower edge of the slope window E (kV/cm)", f"{lo_d:g}"))
            hi = float(_ask("Upper edge of the slope window E (kV/cm)", f"{hi_d:g}"))
            if lo >= hi:
                raise ValueError("the lower edge must be below the upper edge")
            break
        except ValueError as e:
            print(f"  ! {e}, try again")
    args.slope_window = [lo, hi]
    print("Which loops? Condition on E_app, e.g.  20  |  10:20  |  >=15  |  10:15, 20  |  all")
    print(f"  (only loops that reach both {lo:g} and {hi:g} kV/cm give a full-window slope)")
    while True:
        spec = _ask("E_app condition", args.slope_eapp or "all")
        try:
            args.eapp_test = parse_eapp_condition(spec)
            args.slope_eapp = spec
            break
        except ValueError as e:
            print(f"  ! {e}")
    print()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder", type=Path,
                    help="folder with .dat files, or a single .dat file")
    ap.add_argument("output_folder", type=Path)
    ap.add_argument("--strain-channel", default="D1", choices=["D1", "D2", "D3"])
    ap.add_argument("--pol-channel", default="P1", choices=["P1", "P2", "P3"])
    ap.add_argument("--thickness-um", type=float, default=None)
    ap.add_argument("--avg-points", type=int, default=0)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--keep-broken", action="store_true",
                    help="keep measurements with Measurement Status != 0 in the results "
                         "(they are still listed in the broken_measurements CSV)")
    ap.add_argument("--slope-window", nargs=2, type=float, default=[-10.0, 10.0],
                    metavar=("E_LO", "E_HI"),
                    help="field window (kV/cm) where the slope of D [nm] vs V+ [V] is fitted "
                         "(default: -10 10)")
    ap.add_argument("--slope-eapp", default="all", metavar="CONDITION",
                    help="which loops go in the slope table, by applied field X = E_app "
                         "(kV/cm): 20 | 10:20 | >=15 | '10:15, 20' | all (default)")
    ap.add_argument("--interactive", action="store_true",
                    help="list the E_app values in the files and ask for the slope "
                         "window and E_app condition")
    args = ap.parse_args(argv)
    if args.slope_window[0] >= args.slope_window[1]:
        ap.error("--slope-window: E_LO must be smaller than E_HI")
    try:
        args.eapp_test = parse_eapp_condition(args.slope_eapp)
    except ValueError as e:
        ap.error(f"--slope-eapp: {e}")

    if args.input_folder.is_file():
        files = [args.input_folder]
    elif args.input_folder.is_dir():
        files = sorted(p for p in args.input_folder.iterdir() if p.suffix.lower() == ".dat")
    else:
        sys.exit(f"Not found: {args.input_folder}")
    if not files:
        sys.exit(f"No .dat files in {args.input_folder}")
    if args.interactive:
        ask_slope_settings(files, args)
    e_lo, e_hi = args.slope_window
    print(f"Slope window: {e_lo:g} ... {e_hi:g} kV/cm, E_app condition: {args.slope_eapp}")
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
    slopes = allres[args.eapp_test(allres["X_kV_cm"].to_numpy())][SLOPE_COLS]
    slopes.to_csv(o / "slopes_all_files.csv", index=False)
    if not args.no_plots and len(allres):
        plot_s_vs_x_all_files(allres, o / "S_vs_X_all_files.png")
        plot_sr_ratio(allres, fits, "All files", o / "Sr_ratio_all_files.png")
        plot_slopes_all_files(slopes, o / "slopes_all_files.png", e_lo, e_hi, args.slope_eapp)
    n_partial = int((~slopes["slope_window_covered"]).sum())
    print(f"Slopes: {len(slopes)} loop(s) meet the E_app condition"
          + (f", {n_partial} of them do not reach {e_lo:g}/{e_hi:g} kV/cm "
             "(slope_window_covered = False, left out of the plot)" if n_partial else ""))
    try:
        with pd.ExcelWriter(o / "summary_all_files.xlsx") as xw:
            s_vs_x.to_excel(xw, sheet_name="S_vs_X", index=False)
            slopes.to_excel(xw, sheet_name="slopes", index=False)
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
