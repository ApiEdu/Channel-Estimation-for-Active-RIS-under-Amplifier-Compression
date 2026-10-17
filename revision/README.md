# Revision experiments (MDPI Telecom, telecom-4583516)

Python 3 (NumPy, SciPy, Matplotlib). Run from this folder; seeds are fixed in each script.

| Script | Produces | Paper |
|---|---|---|
| `common.py` | amplifier models, Bussgang moments, simulator, bootstrap CI | — |
| `exp_theory.py` | static/fluctuating floor validation (Table 3), K-scaling for rank-one/Rician/Rayleigh, Bayes-MMSE benchmark | Figs. 8, 9 |
| `exp_tracker.py` | tracker uniqueness, accuracy (RMS vs kappa table), finite sample, aging, drift, shape mismatch | Fig. 11 |
| `exp_reversal.py` | optimal training amplification vs direct link, headroom, receiver noise | Fig. 10 |
| `exp_scaling.py` | Fig. 4 (scaling in N and in the Rapp smoothness, with the corrected Eq. 8) | Fig. 4 |
| `exp_tracker_main.py` | Figs. 6–7 (tracker vs fixed calibration; probe ratio) | Figs. 6, 7 |
| `make_figs.py` | the four PNG figures (needs the three JSON files written by the scripts above) | Figs. 8–11 |

Order: `python3 exp_theory.py; python3 exp_tracker.py; python3 exp_reversal.py; python3 make_figs.py; python3 exp_tracker_main.py; python3 exp_scaling.py`
(`exp_tracker_main.py` and `exp_scaling.py` write their own PNGs; `make_figs.py` writes Figs. 8-11.)
(`exp_tracker.py` also writes `e1_curves.npy`). Runtime: about 3 minutes in total.

Figs. 1, 3, 5, 12 are drawn by `plot_matlab_figs.py` from `e_fig1.json`, `e_fig3.json`, `e_fig5.json`, `e_fig12.json`; Fig. 2 by `plot_training_amp.py` from `e_training_amp.json`; Fig. 13 by `plot_spread.py` from `e_spread.json`. The numbers come from the MATLAB scripts in `../matlab` (fig1_nmse_floor, fig2_training_amplitude, fig4_gain_mismatch, fig7_saleh, fig8_element_spread; the JSON files record the run used). Figs. 6-7 are drawn by `plot_tracker_main.py` (called by `exp_tracker_main.py`).

Source of each figure: MATLAB (`../matlab`) for Figs. 1, 2, 3, 5, 12, 13; Python (this folder) for Figs. 4, 6, 7, 8, 9, 10, 11. Exception: `matlab/track_kappa.m` now uses the exact RMS-gain table and the guard window, as in the paper.
`exp_tracker_main.py` regenerates Figs. 6–7 (values in the paper come from it, seed 4242); the MATLAB
`fig5/fig6` scripts give the same behaviour with different random draws, so numbers agree only to within
Monte Carlo error. `matlab/fig1_nmse_floor.m` still uses the previous eps_b in the BLMMSE shrinkage factor
(effect < 0.3 dB on the BLMMSE curve of Fig. 1).
