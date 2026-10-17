"""Fig. 4 regenerated: scaling in N (T=N+1) and in the Rapp smoothness, with the CORRECTED closed form (Eq. 8). seed 707."""
import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C
from common import P, db, cn
rng = np.random.default_rng(707)
p, K = P['p'], P['K']; PT = 10**3.0; s1 = P['sig1_2']
def decomp(s, M=4000, L=400):
    r = np.random.default_rng(99); b = np.sqrt(P['rho_b'])*cn(r, M); v = np.sqrt(s1)*cn(r, M, L)
    u = np.sqrt(PT)*b[:, None]+v; g = C.rapp(u, p, P['Asat'], s); kap = (np.mean(np.conj(u)*g)/(PT*P['rho_b']+s1)).real
    d = g-kap*u; dbar = d.mean(1)
    return kap, np.mean(abs(d)**2), np.mean(abs(dbar)**2), np.mean(abs(d-dbar[:, None])**2)+kap**2*s1
def closed(N, s):
    T = N+1; k, sd, ss, sf = decomp(s); den = k*k*PT*P['rho_b']
    return db((ss+N/T*sf)/den + P['sig2_2']/(K*PT*T*k*k*P['rho_g'])/P['rho_b'])
def sim(N, s, mode, nr=300):
    T = N+1; kap = decomp(s)[0]; err = []; den = []
    for _ in range(nr):
        b = np.sqrt(P['rho_b'])*cn(rng, N); G = C.make_G(rng, K, N, P['rho_g']); hd = np.sqrt(P['rho_d'])*cn(rng, K)
        if mode == 'lin': c = C.observe(b, G, hd, PT, rng, 'linear', None, s, T=T); g = p
        else:
            c = C.observe(b, G, hd, PT, rng, 'rapp', P['Asat'], s, T=T); g = kap if mode == 'cmp' else p
        err.append(np.sum(abs(b-c/g)**2)); den.append(np.sum(abs(b)**2))
    return C.boot_db(err, den, rng, 400)
R = {'N': {}, 'S': {}}
print('== (a) N sweep, T=N+1 ==')
for N in (8, 16, 32, 64):
    r = {m: sim(N, 2, m) for m in ('lin', 'cmp', 'conv')}; r['closed'] = closed(N, 2); R['N'][N] = r
    print(f"N={N:3d} linear {r['lin'][0]:7.2f} comp {r['cmp'][0]:6.2f} [{r['cmp'][1]:.2f},{r['cmp'][2]:.2f}] conv {r['conv'][0]:6.2f} | closed form {r['closed']:6.2f}  diff {r['cmp'][0]-r['closed']:+.2f}")
print('== (b) smoothness sweep, N=16 ==')
for s in (0.5, 1, 2, 4, 8):
    r = {m: sim(16, s, m) for m in ('cmp', 'conv')}; r['closed'] = closed(16, s); R['S'][s] = r
    print(f"varsigma={s:3.1f} comp {r['cmp'][0]:6.2f} [{r['cmp'][1]:.2f},{r['cmp'][2]:.2f}] conv {r['conv'][0]:6.2f} | closed form {r['closed']:6.2f}  diff {r['cmp'][0]-r['closed']:+.2f}")
json.dump(R, open('e_scaling.json', 'w'), default=float, indent=1)
plt.rcParams.update({'font.size': 8, 'axes.grid': True, 'grid.alpha': .3, 'font.family': 'serif', 'mathtext.fontset': 'stix', 'savefig.dpi': 300, 'legend.fontsize': 7})
col = {'lin': '#1f77b4', 'cmp': '#d62728', 'conv': '#9467bd'}
fig, ax = plt.subplots(1, 2, figsize=(6.7, 2.8))
Ns = sorted(R['N']); a = ax[0]
for k, lab, mk in (('lin', 'Linear', 'o'), ('cmp', 'Compensated', '^'), ('conv', 'Conventional', 's')):
    y = [R['N'][n][k][0] for n in Ns]
    a.errorbar(Ns, y, yerr=[[R['N'][n][k][0]-R['N'][n][k][1] for n in Ns], [R['N'][n][k][2]-R['N'][n][k][0] for n in Ns]], fmt='-'+mk, color=col[k], ms=4, capsize=2, label=lab)
a.plot(Ns, [R['N'][n]['closed'] for n in Ns], ':', color='k', lw=1.4, label='Eq. (8)')
a.set_xscale('log', base=2); a.set_xticks(Ns); a.set_xticklabels([str(n) for n in Ns]); a.minorticks_off()
a.set_xlabel('RIS elements $N$ ($T=N+1$)'); a.set_ylabel('NMSE [dB]'); a.set_title('(a) Scaling the surface', loc='left', fontsize=8); a.legend(loc='center left')
a = ax[1]; Ss = sorted(R['S'])
for k, lab, mk in (('cmp', 'Compensated', '^'), ('conv', 'Conventional', 's')):
    y = [R['S'][s][k][0] for s in Ss]
    a.errorbar(Ss, y, yerr=[[R['S'][s][k][0]-R['S'][s][k][1] for s in Ss], [R['S'][s][k][2]-R['S'][s][k][0] for s in Ss]], fmt='-'+mk, color=col[k], ms=4, capsize=2, label=lab)
a.plot(Ss, [R['S'][s]['closed'] for s in Ss], ':', color='k', lw=1.4, label='Eq. (8)')
a.set_xscale('log', base=2); a.set_xticks(Ss); a.set_xticklabels(['0.5', '1', '2', '4', '8']); a.minorticks_off()
a.set_xlabel('Rapp smoothness $\\varsigma$'); a.set_ylabel('NMSE [dB]'); a.set_title('(b) Compression smoothness', loc='left', fontsize=8); a.legend(loc='lower right')
fig.tight_layout(); fig.savefig('../fig_scaling.png'); plt.close(fig)
