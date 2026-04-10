# Bridging continuous and discrete tensor formats of multivariate functions via QTT

[![DOI][doi-badge]][doi-link]
[![arXiv][arxiv-badge]][arxiv-link]
[![MATLAB][code-badge]][code-link]

[doi-badge]: https://img.shields.io/badge/DOI-10.48550/arXiv.XXXX.XXXXX-blue
[doi-link]: https://doi.org/10.48550/arXiv.XXXX.XXXXX

[arxiv-badge]: https://img.shields.io/badge/arXiv-XXXX.XXXXX-red
[arxiv-link]: https://arxiv.org/abs/XXXX.XXXXX

[code-badge]: https://img.shields.io/badge/MATLAB-≥R2022a-blue.svg
[code-link]: https://www.mathworks.com/products/matlab.html

This repository contains the code to reproduce the results of the following paper:

Peter Benner, Boris N. Khoromskij, and Bonan Sun<br>
*Bridging continuous and discrete tensor formats of multivariate functions via QTT*<br>
arXiv preprint, 2026

## Requirements

This code requires:
- **MATLAB** version R2022a or later
- **[Chebfun](https://www.chebfun.org/)** package for Chebyshev approximation
- **[TT-Toolbox](https://github.com/oseledets/TT-Toolbox)** for TT utilities

## Reproducing the Results

0. Download [Chebfun](https://www.chebfun.org/) and [TT-Toolbox](https://github.com/oseledets/TT-Toolbox) and make note of their path locations.

1. Open `setpath.m` and modify the paths at the top to point to your Chebfun and TT-Toolbox installations:
   ```matlab
   path_to_chebfun = '~/Documents/MATLAB/chebfun/';
   path_to_TT_Toolbox = '~/Documents/MATLAB/TT-Toolbox/';
   ```

2. Run:
   ```matlab
   >> reproduce_all            % Fast prototyping (~4 min)
   >> reproduce_all('paper')   % Exact paper results (~20 hours)
   ```

   **Fast prototyping** (default): `d=10, m_max=100, repeat=1, d_values=10:2:20`

   **Paper reproduction** (`'paper'`): `d=20, m_max=300, repeat=20, d_values=10:2:30`

   All parameters can also be customized individually:
   ```matlab
   reproduce_all('mode', 'compute', 'd', 15, 'm_max', 200, 'repeat', 5, 'saveresults', true)
   ```

**CLI / headless mode:** Use `'mode','compute'` to save data without plotting, then `'mode','plot'` to plot from saved data:
```bash
matlab -batch 'reproduce_all("paper", "mode", "compute")'   # headless server
matlab -batch 'reproduce_all("mode", "plot", "saveresults", true)'  # local machine
```

Expected runtime: **~4 minutes** with defaults, **~20 hours** with paper parameters.

## Execution Modes

All scripts (`reproduce_all`, `run_compare_plot`, `run_compare_table`, `run_3d_plots`, `run_scaling`) accept a `'mode'` parameter:

- `'compute&plot'` (default) — compute, save `.mat`, and plot
- `'compute'` — compute and save `.mat` only (no display needed)
- `'plot'` — plot from saved `.mat` data (falls back to `'compute&plot'` if no data exists)

Set `'saveresults', true` to write figures/CSV to disk.

Computed data is saved in `data/`:
- `compare_plot_results.mat` — Experiment 1: single polynomial comparison
- `compare_table_results.mat` — Experiment 2: joint QTT tables
- `3d_plots_results.mat` — Experiment 3: 3D visualizations
- `scaling_results.mat` — Experiment 4: scaling vs. d

## Running Individual Experiments

```matlab
setpath();

% Experiment 1: Single polynomial comparison (Figure 3)
run_compare_plot('repeat', 1, 'saveresults', true, 'd', 20, 'm_max', 300);

% Experiment 2: Joint QTT tables (Tables 2-3)
run_compare_table('repeat', 1, 'saveresults', true, 'd', 20);

% Experiment 3: 3D visualizations (Figure 4)
run_3d_plots('saveresults', true, 'd', 20);

% Experiment 4: Scaling with respect to d (Figure 5)
run_scaling('repeat', 20, 'saveresults', true);
```

## Output

### Figures

All figures are saved to the `figures/` directory:

| File | Description |
|------|-------------|
| [`singlepoly.pdf`](figures/singlepoly.pdf) | Fig. 3: Single polynomial comparison (4 panels × 3 rows: error, rank, runtime) |
| [`biomol.pdf`](figures/biomol.pdf) | Fig. 4(A): Biomolecule potential surface |
| [`biomol_err.pdf`](figures/biomol_err.pdf) | Fig. 4(B): Biomolecule approximation error |
| [`runge.pdf`](figures/runge.pdf) | Fig. 4(C): Runge function surface |
| [`runge_err.pdf`](figures/runge_err.pdf) | Fig. 4(D): Runge approximation error |
| [`wagon.pdf`](figures/wagon.pdf) | Fig. 4(E): Wagon function surface |
| [`wagon_err.pdf`](figures/wagon_err.pdf) | Fig. 4(F): Wagon approximation error |
| [`scaling.pdf`](figures/scaling.pdf) | Fig. 5: Scaling with respect to d (2 columns × 4 rows: error, rank, storage, runtime) |

### Tables

Tables are printed to the MATLAB command window and saved as CSV:

| File | Description |
|------|-------------|
| [`joint_qtt_results.csv`](figures/joint_qtt_results.csv) | Tables 2-3: Joint QTT approximation comparison |

## Note on Figures

Figures in the published paper were generated using TikZ/pgfplots from CSV data. This reproducibility package generates equivalent MATLAB figures directly. The numerical results are identical.

**Note on biomolecule data:** The biomolecule potential data (`biomol.mat`) is stored in native MATLAB format (Chebyshev coefficients and core tensor) and reconstructed as a `chebfun3` object at runtime. This ensures compatibility across different Chebfun versions. The data was generated using [ChebTuck](https://github.com/bonans/ChebTuck). We only provide the processed `biomol.mat` file in the `data/` directory for self-contained and fast reproduction. The original ChebTuck code is available for reference but is not required to run the experiments in this package.

## Code Structure

```
chebqtt/
├── reproduce_all.m          # Main function to run all experiments
├── setpath.m                # Path setup and toolbox verification
├── run_compare_plot.m       # Experiment 1: Single polynomial comparison
├── run_compare_table.m      # Experiment 2: Joint QTT tables
├── run_3d_plots.m           # Experiment 3: 3D visualizations
├── run_scaling.m            # Experiment 4: Scaling with respect to d
├── src/                     # Core source files
│   ├── poly2qtt.m           # Polynomial to QTT conversion
│   ├── ChebTuck2TuckQTT.m   # ChebTuck to QTT-Tucker conversion
│   ├── TuckQTT_eval.m       # QTT-Tucker evaluation
│   └── ...                  # Helper functions
├── data/                    # Input data & computed results
│   ├── biomol.mat           # Biomolecule potential
│   ├── compare_plot_results.mat   # Experiment 1 data (generated)
│   ├── compare_table_results.mat  # Experiment 2 data (generated)
│   ├── 3d_plots_results.mat       # Experiment 3 data (generated)
│   └── scaling_results.mat        # Experiment 4 data (generated)
└── figures/                 # Output figures (generated)
```

## License

MIT License - see [LICENSE.txt](LICENSE.txt)

## Citation

If you use this code, please cite:

```bibtex
@article{chebqtt2025,
  title={Bridging continuous and discrete tensor formats of multivariate functions via QTT},
  author={Benner, Peter and Khoromskij, Boris N. and Sun, Bonan},
  journal={arXiv preprint arXiv:XXXX.XXXXX},
  year={2025}
}
```
