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
python pzm_analysis.py  INPUT_FOLDER  OUTPUT_FOLDER  [options]
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

Every `.dat` file in `INPUT_FOLDER` is processed. A table that fails to parse
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
├── summary_all_files.xlsx         same + one sheet per file
└── <file_name>/
    ├── <file>_results.csv              parameters per loop
    ├── <file>_loops.csv                all loop points (long format: table, t, V, E, D, S, P)
    ├── <file>_instrument_summary.csv   the instrument's own summary table
    ├── <file>_SE.png                   all S‑E loops, coloured by applied field
    ├── <file>_PE.png                   all P‑E loops, coloured by applied field
    ├── <file>_params_vs_E.png          Smax+, Sr, Smax−, 2Pr, Ec vs applied field
    └── <file>_annotated_max.png        largest loop with Smax+, Smax−, Sr marked
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
