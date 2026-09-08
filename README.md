# MATLAB code — Channel Estimation for Active RIS under Amplifier Compression

Reproduces every figure in the paper.

## !! Read this before running !!

**This code has not been executed.** It was written and checked by hand; no
MATLAB or Octave interpreter was available when it was produced. Expect to fix
small errors on the first run. Verify each figure against the numbers in the
paper before you rely on it or release it.

Requires MATLAB R2016b or later (implicit expansion, `pchip` interpolation).
No toolboxes beyond base MATLAB.

## Files

| File | Role |
|---|---|
| `params.m` | All parameters, Table 1 |
| `rapp.m` | Soft-limiting amplifier, eq. (2) |
| `saleh.m` | AM/AM + AM/PM amplifier, eq. (13) |
| `bussgang_coeff.m` | Numerical kappa and sigma_d^2, eq. (4) |
| `sim_nmse.m` | Shared Monte Carlo core for every figure |
| `track_kappa.m` | Algorithm 1, dual-power probing |
| `fig1_nmse_floor.m` … `fig9_scaling.m` | One script per figure |
| `run_all.m` | Runs everything (several hours) |

## Running

Run scripts individually while checking. `run_all.m` only once you trust them.

```matlab
fig1_nmse_floor      % error floor and non-monotonicity
fig5_kappa_tracking  % the proposed algorithm
```

Each script prints its numbers to the console and saves `figN_data.mat`.

## Expected results

| Figure | Key number |
|---|---|
| 1 | knee at IBO = 5 dB; conventional min -13.4 dB then degrades to -0.44 dB |
| 2 | true optimum a = 3.5; penalty 3.4 dB at a = 10 |
| 3 | linear gains 16.9 dB from K=1 to 64; compensated only 2.8 dB |
| 4 | crossover at +3 dB mismatch; 14.1 dB penalty at +10 dB |
| 5 | 10.3 dB gain over fixed calibration at 6 dB drift |
| 6 | probe optimum at 10 dB |
| 7 | phase saturates near 25 deg; 4.0 dB from phase correction |
| 8 | per-element oracle flat; single-gain loses 3.3 dB at sigma_A = 4 dB |
| 9 | compensated floor converges to gamma_d = -6.94 dB, independent of N |

If a figure disagrees with these, the discrepancy is real and must be resolved
before submission — not explained away.

## Two design points that cost dB if you change them

1. **DFT columns 1..N only.** The DC column is reserved for the direct link.
   Including it in the RIS block makes the two collinear and LS becomes
   singular.
2. **Line-of-sight backward link.** A Rayleigh RIS-to-BS link produces
   near-null columns that destroy the conditioning of the LS problem.

## Random seeds

Each script calls `rng` with a fixed seed. Results will still differ slightly
from the paper because the reference implementation was in Python with a
different generator. The trends and the reported figures should match; exact
decimals will not.

---

## Figure export (MDPI Telecom)

`mdpi_style.m` applies the journal's figure requirements and writes the PNG.
Every `figN_*.m` script calls it on its last line, so simply running a script
regenerates the exact file the manuscript includes.

    mdpi_style('fig_nmse_floor_ibo')      % single column, 12 cm
    mdpi_style('fig_scaling', 17)         % full width, two panels

What it sets:

| Item | Value | Why |
|---|---|---|
| Resolution | 300 dpi | MDPI minimum for colour figures |
| Width | 12 cm (17 cm for Fig. 9) | Telecom is single-column |
| Axis font | 11 pt, labels 12 pt | the IEEE two-column sizes read small here |
| Line width | 1.8 pt | thin lines disappear at print size |
| Background | white | avoids grey boxes in the PDF |

**Do not rename the output files.** `paper_telecom.tex` includes them by name;
renaming breaks the build.

To restyle, edit `mdpi_style.m` once rather than the nine scripts.
`exportgraphics` needs MATLAB R2020a or newer; older releases fall back to
`print` automatically.

After regenerating, copy the PNGs into the submission folder and rebuild from
scratch — delete `.aux` and `.bbl` first, or LaTeX will reuse the old figures.
