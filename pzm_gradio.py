#!/usr/bin/env python3
"""
pzm_gradio.py - Web interface for the D-V slope analysis in pzm_analysis.py.

    python pzm_gradio.py                       # opens http://127.0.0.1:7860
    python pzm_gradio.py --data carpeta_dat --out carpeta_resultados --port 7861

Tab "Per file": pick one .dat file (from the data folder or uploaded), set the
field window (kV/cm) and the E_app condition, and look at the slope of
D [nm] vs V+ [V] for every loop, with a check plot of any loop. "Save" writes
<file>_slopes.csv and <file>_slope_check.png for that file to the output
folder, so each file can be saved with its own window.

Tab "Batch": runs the full pzm_analysis.py pipeline (all outputs, all files)
with one window and condition, exactly as the command line does.

The analysis itself lives in pzm_analysis.py. This file only adds the interface.
"""

from __future__ import annotations

import argparse
import contextlib
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gradio as gr

import pzm_analysis as pzm

# compact columns for the on-screen table
VIEW_COLS = ["step", "table", "ramp", "X_kV_cm", "slope_Vlo_V", "slope_Vhi_V",
             "slope_window_covered", "slope_up_nm_per_V", "slope_down_nm_per_V",
             "slope_mean_nm_per_V", "slope_up_R2", "slope_down_R2",
             "slope_up_n", "slope_down_n", "is_last_at_amplitude"]


def list_dat(folder):
    p = Path(folder or ".")
    return sorted(f.name for f in p.iterdir() if f.suffix.lower() == ".dat") if p.is_dir() else []


# --------------------------------------------------------------------------
# Per-file tab
# --------------------------------------------------------------------------

def compute(data_dir, fname, upload, channel, e_lo, e_hi, condition, keep_broken):
    """Analyse one file. Returns everything the tab displays + the state."""
    plt.close("all")        # earlier figures were already rendered by Gradio
    path = Path(upload) if upload else (Path(data_dir) / fname if fname else None)
    if path is None or not path.is_file():
        raise gr.Error("Choose a .dat file from the list or upload one.")
    if e_lo is None or e_hi is None or e_lo >= e_hi:
        raise gr.Error("The lower window edge must be below the upper one.")
    try:
        test = pzm.parse_eapp_condition(condition)
    except ValueError as e:
        raise gr.Error(str(e))

    warnings = []
    out = pzm.analyse_file(path, channel, "P1", None, 0, (e_lo, e_hi),
                           keep_broken, log=warnings.append)
    if out is None:
        raise gr.Error(f"No measurement tables found in {path.name}.")
    res, good, broken, loops, _ = out
    sel = good[test(good["X_kV_cm"].to_numpy())]
    usable = sel[sel["slope_window_covered"]]

    th = good["thickness_um"].iloc[0] if len(good) else float("nan")
    lines = [f"**{path.name}** — {len(res)} loops, {len(broken)} broken "
             f"({'kept' if keep_broken else 'excluded'}), thickness {th:.0f} µm",
             f"Window {e_lo:g} … {e_hi:g} kV/cm = {e_lo * th / 10:.0f} … {e_hi * th / 10:.0f} V; "
             f"fit of **{channel} [nm] vs V+ [V]** on each half-cycle.",
             f"E_app condition `{condition or 'all'}`: {len(sel)} loop(s), "
             f"{len(usable)} reach both window edges."]
    if len(usable):
        last = usable[usable["is_last_at_amplitude"]]
        lines.append(f"Mean slope (last loop at each amplitude): "
                     f"**{last['slope_mean_nm_per_V'].mean():.4f} nm/V** "
                     f"(range {last['slope_mean_nm_per_V'].min():.4f} – "
                     f"{last['slope_mean_nm_per_V'].max():.4f}).")
    lines += [f"⚠ {w.strip()}" for w in warnings]

    fig = pzm.plot_slopes_all_files(sel, None, e_lo, e_hi, condition)
    # tables for the check-plot selector; default = largest usable loop
    pick_from = usable if len(usable) else sel
    choices = [f"table {int(r.table)} - E_app {r.X_kV_cm:.1f} kV/cm"
               for r in pick_from.itertuples()]
    default = None
    if len(pick_from):
        t_def = int(pick_from.loc[pick_from["amplitude_V"].idxmax(), "table"])
        default = next(c for c in choices if c.startswith(f"table {t_def} "))
    state = {"path": str(path), "stem": path.stem, "channel": channel,
             "sel": sel, "loops": loops}
    check = draw_check(state, default)
    view = sel[VIEW_COLS].round(5)
    return ("\n\n".join(lines), view, fig, gr.update(choices=choices, value=default),
            check, state, None)


def draw_check(state, choice):
    if not state or not choice:
        return None
    t = int(choice.split()[1])
    for n in plt.get_fignums()[:-1]:   # keep only the newest (the overview)
        plt.close(n)
    row = state["sel"][state["sel"]["table"] == t].iloc[0]
    loop = state["loops"][state["loops"]["table"] == t]
    return pzm.plot_slope_check(state["stem"], row, loop, None, d_label=state["channel"])


def save_file(state, out_dir, choice):
    """Write <stem>_slopes.csv (+ check plot of the selected loop) for this file."""
    if not state:
        raise gr.Error("Compute first.")
    outdir = Path(out_dir) / state["stem"]
    outdir.mkdir(parents=True, exist_ok=True)
    csv = outdir / f"{state['stem']}_slopes.csv"
    state["sel"][pzm.SLOPE_COLS].to_csv(csv, index=False)
    msg = f"Saved {csv}"
    if choice:
        fig = draw_check(state, choice)
        png = outdir / f"{state['stem']}_slope_check.png"
        fig.savefig(png, dpi=200); plt.close(fig)
        msg += f" and {png.name}"
    return msg, str(csv)


# --------------------------------------------------------------------------
# Batch tab
# --------------------------------------------------------------------------

def run_batch(data_dir, out_dir, channel, e_lo, e_hi, condition, keep_broken, no_plots):
    if e_lo is None or e_hi is None or e_lo >= e_hi:
        raise gr.Error("The lower window edge must be below the upper one.")
    argv = [str(data_dir), str(out_dir), "--strain-channel", channel,
            "--slope-window", str(e_lo), str(e_hi), "--slope-eapp", condition or "all"]
    argv += ["--keep-broken"] if keep_broken else []
    argv += ["--no-plots"] if no_plots else []
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            pzm.main(argv)
    except SystemExit as e:          # argparse errors / no files
        buf.write(f"\nStopped: {e}")
    png = Path(out_dir) / "slopes_all_files.png"
    csv = Path(out_dir) / "slopes_all_files.csv"
    return (buf.getvalue(), str(png) if png.is_file() and not no_plots else None,
            str(csv) if csv.is_file() else None)


# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------

def build(data_dir="carpeta_dat", out_dir="carpeta_resultados"):
    files = list_dat(data_dir)
    with gr.Blocks(title="PZM slope analysis") as app:
        gr.Markdown("## PZM loops: slope of D [nm] vs V+ [V] inside a field window\n"
                    "The window is given in field (kV/cm) and converted to volts with each "
                    "loop's thickness. The slope is fitted separately on the decreasing (V+max → V−max) "
                    "and increasing (V−max → V+max) half-cycles; *mean* is their average.")
        with gr.Row():
            data_in = gr.Textbox(data_dir, label="Data folder (.dat files)")
            out_in = gr.Textbox(out_dir, label="Output folder")
        with gr.Row():
            e_lo = gr.Number(-10, label="Window lower edge E (kV/cm)")
            e_hi = gr.Number(10, label="Window upper edge E (kV/cm)")
            cond = gr.Textbox("all", label="E_app condition",
                              info="20 | 10:20 | >=15 | 10:15, 20 | all  (kV/cm, ±0.05)")
            chan = gr.Dropdown(["D1", "D2", "D3"], value="D1", label="Displacement channel")
            keep = gr.Checkbox(False, label="Keep broken (status ≠ 0)")

        with gr.Tab("Per file"):
            with gr.Row():
                fsel = gr.Dropdown(files, value=files[0] if files else None,
                                   label="File in data folder", scale=3)
                refresh = gr.Button("↻ Refresh list", scale=1)
                upl = gr.File(label="…or upload a .dat (overrides the list)",
                              file_types=[".dat"], type="filepath", scale=3)
            go = gr.Button("Compute slopes", variant="primary")
            info = gr.Markdown()
            with gr.Row():
                overview = gr.Plot(label="Mean slope vs E_app")
                with gr.Column():
                    pick = gr.Dropdown([], label="Loop for the check plot")
                    check = gr.Plot(label="D-V loop with the fitted lines")
            table = gr.Dataframe(label="Slopes (loops meeting the E_app condition)",
                                 interactive=False, wrap=False)
            with gr.Row():
                save = gr.Button("Save this file's slopes to the output folder")
                saved_msg = gr.Markdown()
            saved_csv = gr.File(label="Download slopes CSV")
            state = gr.State(None)

            refresh.click(lambda d: gr.update(choices=list_dat(d)), data_in, fsel)
            go.click(compute, [data_in, fsel, upl, chan, e_lo, e_hi, cond, keep],
                     [info, table, overview, pick, check, state, saved_csv])
            pick.change(draw_check, [state, pick], check)
            save.click(save_file, [state, out_in, pick], [saved_msg, saved_csv])

        with gr.Tab("Batch (all files)"):
            gr.Markdown("Runs the complete `pzm_analysis.py` pipeline on every .dat file in "
                        "the data folder with the settings above, writing all outputs "
                        "(results, plots, Excel) to the output folder — same as the CLI.")
            nop = gr.Checkbox(False, label="No plots (CSV / Excel only, faster)")
            run = gr.Button("Run batch", variant="primary")
            log = gr.Textbox(label="Log", lines=14)
            with gr.Row():
                bplot = gr.Image(label="slopes_all_files.png", type="filepath")
                bcsv = gr.File(label="slopes_all_files.csv")
            run.click(run_batch, [data_in, out_in, chan, e_lo, e_hi, cond, keep, nop],
                      [log, bplot, bcsv])
    return app


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="carpeta_dat", help="folder with .dat files")
    ap.add_argument("--out", default="carpeta_resultados", help="output folder")
    ap.add_argument("--port", type=int, default=7860)
    ap.add_argument("--share", action="store_true", help="public gradio.live link")
    a = ap.parse_args()
    build(a.data, a.out).launch(server_port=a.port, share=a.share)


if __name__ == "__main__":
    main()
