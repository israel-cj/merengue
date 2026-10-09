#!/usr/bin/env python3
"""
pzm_analysis.py - Batch extraction of S-E / P-E loop parameters from aixACCT
TF Analyzer "PiezoResult" (.dat) files (PZM module).

For every measurement table (= one sine period at a given amplitude) it extracts:

  Strain loop (S-E, main figure)
    Smax+  strain at the positive field maximum  (V = V+max)
    Smax-  strain at the negative field maximum  (V = V-max)
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

Outputs (in OUTPUT_FOLDER)
-------
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
    Smax_m, Emax_m = _value_at_extreme(E, S, imin, avg)
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
        "instr_status": meta.get("Measurement Status", np.nan),
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
    # when the same amplitude is measured several times in a row, only the last
    # one (the most stabilised loop) gets is_last_at_amplitude = True
    nxt = np.r_[amp[1:], np.nan]
    res["is_last_at_amplitude"] = amp != nxt
    return res


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
           (r["Emax-_kV_cm"], r["Smax-_%"], "S$_{max}^-$"),
           (0.0, r["Sr_%"], "S$_r$")]
    for x, y, lab in pts:
        ax.plot(x, y, "o", color="crimson", zorder=5)
        ax.annotate(f"{lab} = {y:.4f} %", (x, y), textcoords="offset points",
                    xytext=(8, 8), fontsize=9, color="crimson")
    ax.set_xlabel("E (kV/cm)"); ax.set_ylabel("Strain (%)")
    ax.set_title(f"{stem} - table {int(r['table'])}, E$_{{app}}$ = {r['Eapp_kV_cm']:.1f} kV/cm", fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / f"{stem}_annotated_max.png", dpi=200); plt.close(fig)


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
    loops = pd.concat(loops, ignore_index=True)

    outdir = outroot / stem
    outdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(outdir / f"{stem}_results.csv", index=False)
    loops.to_csv(outdir / f"{stem}_loops.csv", index=False)
    if summary is not None:
        summary.to_csv(outdir / f"{stem}_instrument_summary.csv", index=False)
    if not args.no_plots:
        make_plots(stem, res, loops, outdir)
    print(f"  {path.name}: {len(res)} loops, thickness {res['thickness_um'].iloc[0]:.0f} um, "
          f"Eapp {res['Eapp_kV_cm'].min():.1f}-{res['Eapp_kV_cm'].max():.1f} kV/cm")
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_folder", type=Path)
    ap.add_argument("output_folder", type=Path)
    ap.add_argument("--strain-channel", default="D1", choices=["D1", "D2", "D3"])
    ap.add_argument("--pol-channel", default="P1", choices=["P1", "P2", "P3"])
    ap.add_argument("--thickness-um", type=float, default=None)
    ap.add_argument("--avg-points", type=int, default=0)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args(argv)

    files = sorted(p for p in args.input_folder.iterdir() if p.suffix.lower() == ".dat")
    if not files:
        sys.exit(f"No .dat files in {args.input_folder}")
    args.output_folder.mkdir(parents=True, exist_ok=True)
    print(f"Processing {len(files)} file(s) -> {args.output_folder}")

    allres = [r for f in files if (r := process_file(f, args.output_folder, args)) is not None]
    allres = pd.concat(allres, ignore_index=True)
    allres.to_csv(args.output_folder / "summary_all_files.csv", index=False)
    try:
        with pd.ExcelWriter(args.output_folder / "summary_all_files.xlsx") as xw:
            allres.to_excel(xw, sheet_name="all", index=False)
            for fname, g in allres.groupby("file"):
                g.to_excel(xw, sheet_name=Path(fname).stem[-31:], index=False)
    except ImportError:
        print("  (openpyxl not installed: Excel summary skipped, CSV written)")
    print("Done.")


if __name__ == "__main__":
    main()
