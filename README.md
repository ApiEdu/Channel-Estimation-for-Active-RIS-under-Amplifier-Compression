# Channel Estimation for Active RIS under Amplifier Compression

Code that produces the numbers and figures of the manuscript *Channel Estimation for Active
Reconfigurable Intelligent Surfaces under Amplifier Compression* (submitted to *Telecom*, MDPI).

* MATLAB code (this folder): amplifier models, Monte Carlo core, tracker, and the scripts for Figs. 1, 2, 3, 5, 6, 7, 12, 13.
* `revision/`: Python scripts for the corrected theory, the Bayes benchmark, the tracker robustness study and the
  reversal map (Figs. 4, 6, 7, 8, 9, 10, 11), and the plotting scripts that draw Figs. 1, 2, 3, 5, 12, 13 from the saved numbers.
* `figures/`: the figure files used in the manuscript.

## Requirements
MATLAB R2016b or later (base MATLAB only) or GNU Octave. The scripts were run in MATLAB by the author and
in GNU Octave 8.4. Python 3 with numpy, scipy and matplotlib for `revision/`.

## MATLAB files
| File | Role |
|---|---|
| `params.m` | Simulation parameters (Table 1) |
| `rapp.m`, `saleh.m` | Rapp and Saleh amplifier models (`Asat` may be a scalar or a 1xN vector) |
| `bussgang_coeff.m` | Numerical Bussgang gain and distortion power |
| `sim_nmse.m` | Shared Monte Carlo core |
| `track_kappa.m` | Dual-power gain tracking (Algorithm 1) |
| `fig1_nmse_floor.m` ... `fig9_scaling.m` | One script per figure (see the table below) |
| `mdpi_style.m` | Figure export (300 dpi) |

Run the scripts one at a time; each prints its numbers. `run_all.m` takes several hours.

## Script -> figure in the manuscript
| Script | Manuscript figure | Expected key numbers (Monte Carlo error up to about 0.2 dB) |
|---|---|---|
| `fig1_nmse_floor.m` | Fig. 1 | Conventional LS minimum -13.4 dB at 15 dBm (IBO = 5 dB), -0.44 dB at 45 dBm; compensated floor about -5.5 dB |
| `fig2_training_amplitude.m` | Fig. 2 | True optimum a = 3.5; 3.4 dB penalty at a = 10 |
| `fig3_antenna_scaling.m` | Fig. 3 | Linear model gains 16.7 dB from K = 1 to 64; compensated 2.9 dB |
| `fig4_gain_mismatch.m` | Fig. 5 | Compensation with nominal calibration is worse than none beyond +3 dB mismatch; 14.3 dB penalty at +10 dB |
| `fig5_kappa_tracking.m` | Fig. 6 (MATLAB reproduction) | Estimated gain close to the true gain |
| `fig6_probe_design.m` | Fig. 7 (MATLAB reproduction) | Best probe ratio near 10 dB |
| `fig7_saleh.m` | Fig. 12 | Phase rotation saturates near 25 degrees; phase correction worth 3.8 dB at IBO = 0 dB |
| `fig8_element_spread.m` | Fig. 13 | Per-element oracle flat within 0.2 dB; single nominal gain about 3.2 dB worse at sigma_A = 4 dB |
| `fig9_scaling.m` | superseded | Fig. 4 is now produced by `revision/exp_scaling.py` |

Random-number generators differ between MATLAB, Octave and Python, so independent runs agree to within
Monte Carlo error but not to the last decimal. If a figure disagrees by more than that, please open an issue.

## Notes
* **Vector `Asat`.** `rapp.m` and `saleh.m` divide element-wise (`./`) so that a per-element saturation level works;
  an earlier version used `/` and gave wrong results for per-element spreads (Fig. 13 only).
* **DFT training columns 1..N only.** The DC column is reserved for the direct link; including it makes LS singular.
* **Line-of-sight backward link** is the default; Rician and Rayleigh backward channels are used only in `revision/`.

## Python (`revision/`)
Order: `python3 exp_theory.py; python3 exp_tracker.py; python3 exp_reversal.py; python3 make_figs.py; python3 exp_tracker_main.py; python3 exp_scaling.py`
(about 3 minutes in total). See `revision/README.md` for the source of each figure.
