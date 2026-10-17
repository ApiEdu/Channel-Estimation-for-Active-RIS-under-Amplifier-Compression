import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import common as C
from common import db
plt.rcParams.update({'font.size': 8, 'axes.grid': True, 'grid.alpha': .3, 'legend.fontsize': 7, 'lines.linewidth': 1.4,
                     'font.family': 'serif', 'mathtext.fontset': 'stix', 'savefig.dpi': 300})
th = json.load(open('e_theory.json')); tr = json.load(open('e_tracker.json')); rv = json.load(open('e_reversal.json'))
W = 6.7
col = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e']

# ---------------- Fig A: corrected floor + K scaling -------------------------
fig, ax = plt.subplots(1, 2, figsize=(W, 3.3))
for s1db, c, mk in ((-70, col[0], 'o'), (-30, col[1], 's')):
    rows = [r for r in th['val'] if r['geom'] == 'rank1' and r['s1_dBm'] == s1db and r['PT_dBm'] == 30 and r['K'] == 16]
    T = [r['T'] for r in rows]
    ax[0].errorbar(T, [r['sim'][0] for r in rows], yerr=[[r['sim'][0]-r['sim'][1] for r in rows], [r['sim'][2]-r['sim'][0] for r in rows]],
                   fmt=mk, color=c, ms=4, capsize=2, label=f'Simulation, $\\sigma_1^2={s1db}$ dBm')
    ax[0].plot(T, [r['new'] for r in rows], '-', color=c)
    ax[0].plot(T, [r['old'] for r in rows], '--', color=c, alpha=.8)
ax[0].plot([], [], 'k-', label='Eq. (8), corrected'); ax[0].plot([], [], 'k--', label='Previous form $(N/T)\\gamma_d$')
ax[0].set_xscale('log'); ax[0].minorticks_off(); ax[0].set_xticks([17, 34, 68, 136]); ax[0].set_xticklabels(['17', '34', '68', '136'])
ax[0].set_xlabel('Training slots $T$'); ax[0].set_ylabel('Compensated NMSE [dB]'); ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.27), ncol=2, frameon=False, fontsize=6); ax[0].set_title('(a)', loc='left')
Ks = [2, 4, 8, 16, 32, 64, 128]
for (lab, mk, c), key in zip((('Rank-one LoS', 'o', col[0]), ('Rician 10 dB', 's', col[2]), ('Rician 0 dB', '^', col[3]), ('Rayleigh', 'd', col[1])),
                             ('rank-one LoS', 'Rician 10 dB', 'Rician 0 dB', 'Rayleigh')):
    rows = th['K'][key]
    ax[1].errorbar(Ks, [r['sim'][0] for r in rows], yerr=[[r['sim'][0]-r['sim'][1] for r in rows], [r['sim'][2]-r['sim'][0] for r in rows]],
                   fmt=mk, color=c, ms=4, capsize=1.5, label=lab)
    if rows[0]['th'] is not None and key in ('rank-one LoS', 'Rayleigh'): ax[1].plot(Ks, [r['th'] for r in rows], '-', color=c, lw=.9)
ax[1].plot(Ks, [l[0] for l in th['K']['linear']], 'k--', label='Linear model')
ax[1].set_xscale('log', base=2); ax[1].set_xlabel('BS antennas $K$'); ax[1].set_ylabel('Compensated NMSE [dB]'); ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.27), ncol=3, frameon=False, fontsize=6); ax[1].set_title('(b)', loc='left')
fig.tight_layout(); fig.savefig('../fig_floor_decomp.png'); plt.close(fig)

# ---------------- Fig B: Bayes benchmark ------------------------------------
fig, ax = plt.subplots(figsize=(3.4, 3.3)); B = th['bayes']; x = [r['PT_dBm'] for r in B]
for key, lab, mk, c in (('conv', 'Conventional LS', 's', col[1]), ('comp', 'Compensated LS', '^', col[2]), ('blmmse', 'Bussgang-LMMSE', 'o', col[0]), ('bayes', 'Bayes-MMSE, true Rapp', 'x', 'k')):
    ax.errorbar(x, [r[key][0] for r in B], yerr=[[r[key][0]-r[key][1] for r in B], [r[key][2]-r[key][0] for r in B]], fmt='-'+mk, color=c, ms=3.5, capsize=1.5, lw=1, label=lab)
ax.set_xlabel('Pilot power $P_T$ [dBm]'); ax.set_ylabel('NMSE [dB]'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=6.5)
fig.tight_layout(); fig.savefig('../fig_bayes_benchmark.png'); plt.close(fig)

# ---------------- Fig C: tracker robustness ---------------------------------
x_, kap_, rms_ = np.load('e1_curves.npy', allow_pickle=True)[0].real, np.load('e1_curves.npy')[1].real, np.load('e1_curves.npy')[2].real
rho = 10.0; lg = np.log(x_)
def curve(g):
    gi = lambda v: np.exp(np.interp(np.log(v), lg, np.log(g))); ib = x_[x_*rho < x_[-1]]; return ib, gi(ib)/gi(ib*rho)
fig, ax = plt.subplots(2, 3, figsize=(W, 4.4)); ax = ax.ravel()
for g, lab, c, ls in ((rms_, 'RMS gain', col[0], '-'), (kap_, '$|\\kappa|/p$ (previous)', col[1], '--')):
    ib, r = curve(g); ax[0].plot(db(ib), r, ls, color=c, label=lab)
ax[0].axvspan(-12, 8, color='gray', alpha=.12, label='Guard window'); ax[0].set_xlim(-30, 30)
ax[0].set_xlabel('Back-off $\\mathrm{IBO}_2$ [dB]'); ax[0].set_ylabel('Probe ratio $r$'); ax[0].legend(fontsize=6, loc='upper left'); ax[0].set_title('(a)', loc='left')
for kind, lab, c, mk in (('rms', 'RMS table', col[0], 'o'), ('kappa', 'Kappa table', col[1], 's')):
    rows = tr['acc_'+kind]; ib = [r['ibo_true_dB'] for r in rows]
    ax[1].plot(ib, [r['penalty_dB'] for r in rows], '-'+mk, color=c, ms=3, label=lab)
ax[1].axhline(0, color='k', lw=.6); ax[1].set_xlabel('True back-off [dB]'); ax[1].set_ylabel('NMSE vs oracle $\\kappa$ [dB]'); ax[1].legend(fontsize=6); ax[1].set_title('(b)', loc='left')
fs = tr['fs']; Ns = sorted(int(k) for k in fs)
ax[2].plot(Ns, [fs[str(n)]['ibo_std_dB'] for n in Ns], '-o', color=col[0], ms=3, label='std of $\\widehat{\\mathrm{IBO}}$ [dB]')
ax[2].plot(Ns, [100*fs[str(n)]['kap_relstd'] for n in Ns], '-s', color=col[1], ms=3, label='std of $\\hat\\kappa/\\kappa$ [%]')
ax[2].set_xscale('log', base=2); ax[2].set_xlabel('RIS elements $N$'); ax[2].legend(fontsize=6); ax[2].set_title('(c)', loc='left')
ag = tr['age']; ra = sorted(ag, key=float, reverse=True)
ax[3].plot([float(k) for k in ra], [ag[k]['penalty_dB'] for k in ra], '-o', color=col[0], ms=3); ax[3].invert_xaxis()
ax[3].set_xlabel('Channel correlation $\\rho$ between probes'); ax[3].set_ylabel('Loss vs oracle [dB]'); ax[3].set_title('(d)', loc='left')
dr = tr['drift']; ad = tr['adrift']
ax[4].plot([float(k) for k in sorted(dr, key=float)], [dr[k]['penalty_dB'] for k in sorted(dr, key=float)], '-o', color=col[0], ms=3, label='Channel-gain drift')
ax[4].plot([float(k) for k in sorted(ad, key=float)], [ad[k]['penalty_dB'] for k in sorted(ad, key=float)], '-s', color=col[2], ms=3, label='$A_{sat}$ drift')
ax[4].axhline(0, color='k', lw=.6); ax[4].set_xlabel('Drift between probes [dB]'); ax[4].set_ylabel('Loss vs oracle [dB]'); ax[4].legend(fontsize=6); ax[4].set_title('(e)', loc='left')
sh = tr['shape']; ss = [1, 2, 3, 5, 10]
ax[5].plot(ss, [sh[str(s)]['penalty_dB'] for s in ss], '-o', color=col[0], ms=3); ax[5].set_xscale('log'); ax[5].set_xticks(ss); ax[5].set_xticklabels([str(s) for s in ss])
ax[5].axhline(0, color='k', lw=.6); ax[5].set_xlabel('True smoothness $\\varsigma$ (table: 2)'); ax[5].set_ylabel('Loss vs oracle [dB]'); ax[5].set_title('(f)', loc='left')
fig.tight_layout(); fig.savefig('../fig_tracker_robust.png'); plt.close(fig)

# ---------------- Fig D: reversal map ---------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(W, 2.9))
specs = (('direct_link', lambda x: 10*np.log10(x/1e-6), 'Direct link $\\varrho_d/\\varrho_b$ [dB]'),
         ('headroom', lambda x: 20*np.log10(x/0.10), '$A_{sat}$ relative to nominal [dB]'),
         ('rx_noise', lambda x: 10*np.log10(x), 'Receiver noise $\\sigma_2^2$ [dBm]'))
for a, (k, f, lab), t in zip(ax, specs, 'abc'):
    xs = [f(r['x']) for r in rv[k]]
    a.plot(xs, [r['a_star'] for r in rv[k]], '-o', color=col[0], ms=3.5, label='Optimal $a^*$'); a.axhline(10, color='k', ls='--', lw=.8)
    a.set_ylim(0, 11); a.set_xlabel(lab); a.set_title(f'({t})', loc='left')
    a2 = a.twinx(); a2.plot(xs, [r['pen_amax'] for r in rv[k]], '-s', color=col[1], ms=3, label='Loss at $a=10$'); a2.grid(False)
    a.set_ylabel('Optimal $a^*$', color=col[0]); a2.set_ylabel('Loss at $a_{max}$ [dB]', color=col[1])
    if k == 'direct_link': hh = a.get_lines()[:1] + a2.get_lines()[:1] + a.get_lines()[1:2]
fig.legend(hh, ['Optimal $a^*$ (left axis)', 'Loss at $a=10$ (right axis)', '$a_{max}=10$'], loc='lower center', ncol=3, frameon=False, fontsize=7)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig('../fig_reversal_map.png'); plt.close(fig)
print('figures written')
