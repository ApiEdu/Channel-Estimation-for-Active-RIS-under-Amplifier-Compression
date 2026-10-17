import json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = json.load(open('e_tracker_main.json'))
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .3, 'font.family': 'serif', 'mathtext.fontset': 'stix', 'savefig.dpi': 300})
devs = sorted(R['dev'], key=float)
fig, ax = plt.subplots(figsize=(4.7, 3.2))
for k, lab, c, mk, ls in (('nom', 'Fixed calibration at nominal $A_{sat}$', '#d62728', 's', '-'), ('trk', 'Proposed dual-power tracking', '#ff7f0e', 'o', '-'), ('orc', 'Oracle $\\kappa$', '#2ca02c', '^', '--')):
    y = [R['dev'][d][k][0] for d in devs]; xd = [float(d) for d in devs]
    lo = [R['dev'][d][k][0]-R['dev'][d][k][1] for d in devs]; hi = [R['dev'][d][k][2]-R['dev'][d][k][0] for d in devs]
    ax.errorbar(xd, y, yerr=[lo, hi], fmt=ls+mk, color=c, ms=4, capsize=2, lw=1.4, label=lab)
ax.set_xlabel('$A_{sat}$ deviation from nominal [dB]'); ax.set_ylabel('NMSE at $P_T=30$ dBm [dB]'); ax.legend(fontsize=7, loc='upper left')
ax.set_xlim(-7, 7); fig.tight_layout(); fig.savefig('../fig_kappa_tracking.png'); plt.close(fig)
vs = sorted(R['probe'], key=float)
fig, ax = plt.subplots(figsize=(4.7, 3.7))
y = [R['probe'][v]['trk'][0] for v in vs]
ax.errorbar([float(v) for v in vs], y, yerr=[[R['probe'][v]['trk'][0]-R['probe'][v]['trk'][1] for v in vs], [R['probe'][v]['trk'][2]-R['probe'][v]['trk'][0] for v in vs]], fmt='-o', color='#1f77b4', ms=4, capsize=2, label='NMSE (left axis)')
ax.set_xlabel('Probe power ratio $\\varrho=P_2/P_1$ [dB]'); ax.set_ylabel('NMSE [dB]')
a2 = ax.twinx(); a2.plot([float(v) for v in vs], [R['probe'][v]['ibo_bias_dB'] for v in vs], '--s', color='#2ca02c', ms=4, label='IBO bias (right axis)'); a2.set_ylabel('Bias of estimated IBO [dB]'); a2.grid(False)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
ax.set_xlim(2, 27)
fig.tight_layout(); fig.savefig('../fig_probe_design.png'); plt.close(fig)
