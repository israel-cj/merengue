# PZM loop analysis

Batch extraction of strain (S‑E) and polarization (P‑E) loop parameters from
aixACCT TF Analyzer **PiezoResult** `.dat` files (PZM module).

Point it at a folder of `.dat` files. It writes one results table for all of them,
plus CSVs and plots for each file, to an output folder.

---

## Installation

```bash
# optional: create a virtual environment
python -m venv .venv
# Windows:  .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate

pip install -r requirements.txt
```

Requires Python 3.8 or newer.

## Usage

```bash
python pzm_analysis.py  INPUT_FOLDER_OR_FILE  OUTPUT_FOLDER  [options]
```

e.g.

```bash
python pzm_analysis.py carpeta_dat carpeta_resultados
```

Example (Windows):

```bash
python pzm_analysis.py "C:\Users\P313565\Py\Hyst\data" "C:\Users\P313565\Py\Hyst\results"
```

| Option | Default | Meaning |
|---|---|---|
| `--strain-channel {D1,D2,D3}` | `D1` | Displacement column used for strain |
| `--pol-channel {P1,P2,P3}` | `P1` | Polarization column used |
| `--thickness-um VALUE` | from file | Override the sample thickness (µm) |
| `--avg-points N` | `0` | Average ±N points around V+max / V−max (noise reduction) |
| `--no-plots` | off | Write only the CSV/Excel files |
| `--keep-broken` | off | Keep measurements with `Measurement Status != 0` |
| `--slope-window E_LO E_HI` | `-10 10` | Field window (kV/cm) where the slope of D [nm] vs V+ [V] is fitted |
| `--slope-eapp CONDITION` | `all` | Which loops go in the slope table, by applied field E_app (see below) |
| `--interactive` | off | List the E_app values found in the files, then ask for the window and the condition |

Every `.dat` file in `INPUT_FOLDER` is processed. A single `.dat` file can be
given instead of a folder. A table that fails to parse
is reported and skipped, and the rest of the run continues.

---

## What is extracted

Each `.dat` file has one summary table, then one block per measurement. Each block
is one full sine period at a given amplitude and has its own metadata (thickness,
area, amplitude…). **Every block is analysed separately**, so the number of
amplitudes and their order (ramp up, ramp down, repeated voltages) can differ
between files.

Unit conversion uses the thickness stored in each block:

* `E [kV/cm] = V / thickness`
* `S [%] = D / thickness × 100`

### Strain loop (S‑E)

| Column | Definition |
|---|---|
| `Smax+_%` | Strain at the positive field maximum (V = V+max) |
| `Smax-_%` | Strain at the negative field maximum (V = V−max) |
| `Sr_%` | Remanent strain: strain where V crosses 0 on the descending branch (after V+max) |
| `S_end_%` | Strain when V returns to 0 at the end of the cycle |
| `Spp_%` | Peak‑to‑peak strain (max − min) |
| `Smin_%`, `E_at_Smin_kV_cm` | Minimum strain in the loop and the field where it occurs |
| `d33*_pm_V` | Large‑signal d33\* = Smax+ / Emax+ |

> **Note:** once the sample switches at high field, the strain *at* V−max can
> be positive, and the real negative dip occurs before V−max. `Smax-_%` follows the
> figure definition (strain at V−max). `Smin_%` gives the true minimum.

### Polarization loop (P‑E)

| Column | Definition |
|---|---|
| `Ps+_uC_cm2`, `Ps-_uC_cm2` | Polarization at V+max / V−max |
| `Pr+_uC_cm2`, `Pr-_uC_cm2` | Polarization at V = 0 (descending / ascending branch) |
| `2Pr_uC_cm2` | Pr+ − Pr− |
| `Ec+_kV_cm`, `Ec-_kV_cm` | Field where P crosses 0 (ascending / descending branch) |

Zero crossings are found by linear interpolation between data points.

### Slope of D [nm] vs V+ [V] in a field window

For every loop, a straight line `D = slope · V + b` is fitted (least squares)
to the raw `V+ [V]` and `D1 [nm]` columns (or the `--strain-channel` chosen),
using only the points with `E_LO ≤ E ≤ E_HI` (default −10 … +10 kV/cm). The
window is given in field and converted to volts with each table's thickness
(e.g. ±10 kV/cm = ±710 V at 710 µm). The fit is done separately on each half
of the loop:

| Column | Definition |
|---|---|
| `slope_down_nm_per_V` | Decreasing field half‑cycle (V+max → V−max) |
| `slope_up_nm_per_V` | Increasing field half‑cycle (V−max → V+max) |
| `slope_mean_nm_per_V` | Mean of the two |
| `slope_Vlo_V`, `slope_Vhi_V` | The window in volts for that table |
| `slope_up_R2`, `slope_down_R2`, `slope_*_n` | Fit quality and number of points used |
| `slope_window_covered` | `False` when the loop does not reach both window edges (E_app smaller than the window). The slope is then fitted over a shorter range. These loops stay in the CSV but are left out of the plot |

1 nm/V = 1000 pm/V. Because the thickness cancels (S/E = D/V), this is the same
number as the small‑window d33 from the S‑E loop.

At high field, after the sample switches, the increasing‑field half can be far
from linear inside the window (low R²). Check `*_slope_check.png`, and use the
`up` / `down` columns separately if the mean is not meaningful.

**Choosing the loops (`--slope-eapp`)**, by applied field E_app in kV/cm
(matched within ±0.05 kV/cm; parts can be combined with commas):

| Condition | Selects |
|---|---|
| `all` | every loop (default) |
| `20` | E_app = 20 |
| `10:20` or `10-20` | 10 ≤ E_app ≤ 20 |
| `>=15`, `>15`, `<=12`, `<12` | comparison |
| `"10:15, 20"` | 10 ≤ E_app ≤ 15 **or** E_app = 20 |

```bash
python pzm_analysis.py carpeta_dat carpeta_resultados --slope-window -10 10 --slope-eapp ">=12"
python pzm_analysis.py carpeta_dat carpeta_resultados --interactive
python pzm_analysis.py carpeta_dat/PNZT-10P4-7_P120C_35KVcm_5to20.dat carpeta_resultados --slope-window -5 5
```

### Web interface (Gradio)

```bash
pip install gradio
python pzm_gradio.py            # then open http://127.0.0.1:7860
python pzm_gradio.py --data carpeta_dat --out carpeta_resultados --port 7861
```

* **Per file** tab: pick a `.dat` from the data folder (or upload one), set the
  window (kV/cm), the E_app condition and the channel, then press *Compute*.
  It shows the slope table, the slope‑vs‑E_app plot, and a D‑V check plot for
  any loop you pick. *Save* writes `<file>_slopes.csv` and
  `<file>_slope_check.png` for that file only, so each file can keep its own window.
* **Batch** tab: runs the whole `pzm_analysis.py` pipeline on every file with
  the settings shown, the same as the command line.

#### Columns of the "Slopes (loops meeting the E_app condition)" table

One row per loop (measurement table) in the file whose E_app meets the
condition. Broken measurements (status ≠ 0) are left out unless *Keep broken*
is ticked. Values are rounded to 5 decimals on screen. The saved CSV has full
precision and a few more columns (file, branch, amplitude, thickness, E at
V±max, window in kV/cm, status).

| Column | Unit | Meaning |
|---|---|---|
| `step` | – | Position of the loop in the measurement sequence of the file (1 = first measured) |
| `table` | – | Table number in the `.dat` file (same number as in the instrument software). Use it to pick the loop for the check plot |
| `ramp` | – | `up` if E_app is higher than in the previous loop, `down` if lower, `repeat` if equal |
| `X_kV_cm` | kV/cm | Applied field E_app = amplitude / thickness. This is what the E_app condition is tested on |
| `slope_Vlo_V`, `slope_Vhi_V` | V | The field window converted to volts for this loop: E_LO × thickness and E_HI × thickness (e.g. −710 / +710 V for ±10 kV/cm at 710 µm). Only points with V+ in this range are fitted |
| `slope_window_covered` | True/False | `True` if the loop reaches both window edges. `False` (E_app smaller than the window) means the fit only used the part of the window that the loop reaches. Treat these slopes with care. They are not drawn in the plot |
| `slope_up_nm_per_V` | nm/V | Slope of the straight line D = slope · V + b fitted on the **increasing‑voltage** half of the loop (V−max → V+max) inside the window |
| `slope_down_nm_per_V` | nm/V | Same, on the **decreasing‑voltage** half (V+max → V−max) |
| `slope_mean_nm_per_V` | nm/V | Average of the up and down slopes. This is the value shown in the slope‑vs‑E_app plot. 1 nm/V = 1000 pm/V |
| `slope_up_R2`, `slope_down_R2` | 0–1 | R² of each fit: how close the loop is to a straight line in the window (1 = perfectly linear). Low values (e.g. < 0.9) mean the half‑cycle is curved there (typical after switching at high field), so its slope and the mean are less meaningful |
| `slope_up_n`, `slope_down_n` | points | Number of data points used in each fit (points of that half‑cycle inside the window) |
| `is_last_at_amplitude` | True/False | `False` when the next loop was measured at the same amplitude. Only `True` rows (the last, most stable loop at each amplitude) are used in the plot and in the mean slope shown above the table |

### Bookkeeping columns

| Column | Meaning |
|---|---|
| `step` | Measurement order in the file |
| `ramp` | `up`, `down` or `repeat` compared with the previous amplitude |
| `is_last_at_amplitude` | `False` if the next measurement uses the same amplitude. Filter on `True` to keep only the last (most stable) loop |
| `instr_*` | Values the instrument computed itself (Dvmax±, Pr+, Vc+, status), for cross‑checking |

---

## Output structure

```
OUTPUT_FOLDER/
├── summary_all_files.csv          one row per loop, all files
├── summary_all_files.xlsx         same + one sheet per file (and a "slopes" sheet)
├── slopes_all_files.csv / .png    D‑V slope (nm/V) in the field window vs E_app
└── <file_name>/
    ├── <file>_results.csv              parameters per loop
    ├── <file>_loops.csv                all loop points (long format: table, t, V, E, D, S, P)
    ├── <file>_instrument_summary.csv   the instrument's own summary table
    ├── <file>_SE.png                   all S‑E loops, coloured by applied field
    ├── <file>_PE.png                   all P‑E loops, coloured by applied field
    ├── <file>_params_vs_E.png          Smax+, Sr, Smax−, 2Pr, Ec vs applied field
    ├── <file>_annotated_max.png        largest loop with Smax+, Smax−, Sr marked
    ├── <file>_slopes.csv               slope per loop meeting the E_app condition
    └── <file>_slope_check.png          largest selected loop with the two fitted lines
```

Use `*_annotated_max.png` to check visually that the points are picked correctly.
Use `*_loops.csv` to remake publication figures (e.g. with plotnine/ggplot)
without reparsing the `.dat` files.

### Replotting example

```python
import pandas as pd
from plotnine import ggplot, aes, geom_path, labs, theme_bw

loops = pd.read_csv("results/PI152_3_sine/PI152_3_sine_loops.csv")
res   = pd.read_csv("results/PI152_3_sine/PI152_3_sine_results.csv")
loops = loops.merge(res[["table", "Eapp_kV_cm"]], on="table")

(ggplot(loops, aes("E_kV_cm", "S_%", group="table", color="Eapp_kV_cm"))
 + geom_path() + labs(x="E (kV/cm)", y="Strain (%)") + theme_bw())
```

---

## Validation

The script was tested on 5 files (PI152_3, PI152_7, PNZT‑N4‑P4‑1/2/3) with
16–41 loops each and thicknesses of 780–1000 µm. The extracted Smax+ and Pr+
agree with the instrument's own `Dvmax+` and `Pr+` to < 0.01 nm and < 0.0001 µC/cm².
Smax− agrees with `Dvmax-` to within ~1 nm.

## Differences from the original notebook

* Thickness is read from each file. The notebook hard‑coded 900 µm (`/90`, `/9000`),
  which is wrong for samples of other thicknesses.
* Loops are not grouped by detecting amplitude drops. Each measurement block is
  kept as its own loop, with `step`, `ramp` and `is_last_at_amplitude` labels.
