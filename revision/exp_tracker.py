"""Tracker robustness: uniqueness, finite-sample, channel aging, shape mismatch.  seed 2026."""
import numpy as np, json, sys
import common as C
from common import P, db, cn
rng = np.random.default_rng(2026)
p, N, K = P['p'], P['N'], P['K']; rho = 10**(P['varrho_dB']/10)
PT2 = 10**3.0; PT1 = PT2/rho
tabs = {}
def tab(model, s): 
    k = (model, s)
    if k not in tabs: tabs[k] = C.table(model, s)
    return tabs[k]

# ---- E1: uniqueness / sensitivity of r(IBO2) for both candidate tables -------
x, kap, rms = tab('rapp', 2)
def ratio_curve(g, rho):
    # r(IBO2) = g(IBO2)/g(rho*IBO2); evaluate on the grid by interpolation in log domain
    lg = np.log(x); gi = lambda v: np.exp(np.interp(np.log(v), lg, np.log(g)))
    ib = x[(x*rho < x[-1])]
    return ib, gi(ib)/gi(ib*rho)
out = {}
for name, g in (('kappa', np.abs(kap)), ('rms', rms)):
    ib, r = ratio_curve(g, rho)
    S = np.gradient(np.log(r), np.log(ib))
    out[name] = dict(monotone=bool(np.all(np.diff(r) > 0)), rmin=float(r.min()), rmax=float(r.max()),
                     S_max=float(S.max()), window_dB=[float(db(ib[S > 0.05].min())), float(db(ib[S > 0.05].max()))])
    print(name, out[name])
np.save('e1_curves.npy', np.array([x, np.abs(kap), rms]))

# ---- helpers -----------------------------------------------------------------
GUARD = (-12.0, 8.0)      # dB: informative window; estimates are clipped into it (safeguard)
def inv_ratio(r, g, guard=True):
    ib, rr = ratio_curve(g, rho)
    v = float(np.exp(np.interp(r, rr, np.log(ib))))           # monotone -> unique root (clipped at table ends)
    if guard: v = float(np.clip(v, 10**(GUARD[0]/10), 10**(GUARD[1]/10)))
    return v

def one_run(A_true, model='rapp', s_true=2, model_tab='rapp', s_tab=2, table_kind='rms', rho_age=1.0,
            drift_dB=0.0, A_drift_dB=0.0, Nn=N, Kk=K, geom='rank1', T=None):
    b1 = np.sqrt(P['rho_b'])*cn(rng, Nn)
    # probe 2 (operating) uses b1; probe 1 (low power) sees an aged / drifted channel
    b2 = rho_age*b1+np.sqrt(1-rho_age**2)*np.sqrt(P['rho_b'])*cn(rng, Nn)
    b2 = b2*10**(drift_dB/20)
    G = C.make_G(rng, Kk, Nn, P['rho_g'], geom); hd = np.sqrt(P['rho_d'])*cn(rng, Kk)
    c_op = C.observe(b1, G, hd, PT2, rng, model, A_true, s_true, T=T)
    c_lo = C.observe(b2, G, hd, PT1, rng, model, A_true*10**(A_drift_dB/20), s_true, T=T)
    r = np.linalg.norm(c_op)/np.linalg.norm(c_lo)
    xt, kt, rt = tab(model_tab, s_tab)
    g = rt if table_kind == 'rms' else np.abs(kt)
    ib_hat = inv_ratio(r, g) if True else None
    # kappa from the table at the estimated IBO
    kap_hat = p*np.interp(np.log(ib_hat), np.log(xt), kt.real)+1j*p*np.interp(np.log(ib_hat), np.log(xt), kt.imag)
    su2 = PT2*P['rho_b']+P['sig1_2']
    kap_true = C.moments(model, su2, p, A_true, s_true)[0]
    ib_true = C.ibo(A_true, p, su2)
    e_hat = np.sum(np.abs(b1-c_op/kap_hat)**2); e_or = np.sum(np.abs(b1-c_op/kap_true)**2)
    return ib_hat, ib_true, kap_hat, kap_true, e_hat, e_or, np.sum(np.abs(b1)**2)

def sweep(A_list, nr=300, **kw):
    rows = []
    for A in A_list:
        res = np.array([one_run(A, **kw) for _ in range(nr)], dtype=object)
        ibh = np.array([r[0] for r in res], float); ibt = res[0][1]
        kh = np.array([abs(r[2]) for r in res]); kt = abs(res[0][3])
        eh = np.array([r[4] for r in res], float); eo = np.array([r[5] for r in res], float); dn = np.array([r[6] for r in res], float)
        n_h = C.boot_db(eh, dn, rng, 400); n_o = C.boot_db(eo, dn, rng, 400)
        rows.append(dict(A=A, ibo_true_dB=float(db(ibt)), ibo_med_dB=float(db(np.median(ibh))),
                         ibo_std_dB=float(np.std(db(ibh))), kap_relstd=float(np.std(kh)/kt), kap_relbias=float(np.mean(kh)/kt-1),
                         nmse_trk=n_h, nmse_orc=n_o, penalty_dB=float(n_h[0]-n_o[0])))
    return rows
A_list = [0.02, 0.03, 0.05, 0.07, 0.10, 0.14, 0.20, 0.30]
R = {}
print('\nE2 accuracy (N=K=16): rms-table vs kappa-table'); 
for kind in ('rms', 'kappa'):
    R['acc_'+kind] = sweep(A_list, table_kind=kind)
    for r in R['acc_'+kind]:
        print(f"{kind:5s} IBO={r['ibo_true_dB']:6.1f} dB est={r['ibo_med_dB']:6.1f} sd={r['ibo_std_dB']:4.2f} dB | kappa relbias {100*r['kap_relbias']:6.2f}% sd {100*r['kap_relstd']:5.2f}% | NMSE trk {r['nmse_trk'][0]:6.2f} [{r['nmse_trk'][1]:.2f},{r['nmse_trk'][2]:.2f}] orc {r['nmse_orc'][0]:6.2f} pen {r['penalty_dB']:.2f}")
print('\nE3 finite-sample: IBO std vs N (A=0.10, IBO=-10 dB)')
R['fs'] = {}
for Nn in (8, 16, 32, 64, 128):
    rr = sweep([0.10], nr=300, Nn=Nn, Kk=16)[0]; R['fs'][Nn] = rr
    print(f"N={Nn:3d} ibo sd {rr['ibo_std_dB']:.2f} dB  kappa relstd {100*rr['kap_relstd']:.2f}%  penalty {rr['penalty_dB']:.2f} dB")
print('\nE4 aging (A=0.10): probe-1 channel correlation rho')
R['age'] = {}
for ra in (1.0, 0.99, 0.95, 0.9, 0.8, 0.6):
    rr = sweep([0.10], nr=300, rho_age=ra)[0]; R['age'][ra] = rr
    print(f"rho={ra:4.2f} IBO est {rr['ibo_med_dB']:6.2f} (true {rr['ibo_true_dB']:.2f}) sd {rr['ibo_std_dB']:.2f}  kappa bias {100*rr['kap_relbias']:.2f}% sd {100*rr['kap_relstd']:.2f}% penalty {rr['penalty_dB']:.2f} dB")
print('\nE4b gain drift between probes (large-scale power drift, dB)')
R['drift'] = {}
for dd in (-3, -1, 0, 1, 3):
    rr = sweep([0.10], nr=300, drift_dB=dd)[0]; R['drift'][dd] = rr
    print(f"drift={dd:+d} dB IBO est {rr['ibo_med_dB']:6.2f} kappa bias {100*rr['kap_relbias']:.2f}% penalty {rr['penalty_dB']:.2f} dB")
print('\nE4c Asat drift between probes')
R['adrift'] = {}
for dd in (-2, -1, 0, 1, 2):
    rr = sweep([0.10], nr=300, A_drift_dB=dd)[0]; R['adrift'][dd] = rr
    print(f"Asat drift={dd:+d} dB kappa bias {100*rr['kap_relbias']:.2f}% penalty {rr['penalty_dB']:.2f} dB")
print('\nE5 shape mismatch: table varsigma=2 (Rapp), true varsigma varies, A=0.10')
R['shape'] = {}
for st in (1, 2, 3, 5, 10):
    rr = sweep([0.10], nr=300, s_true=st, s_tab=2)[0]; R['shape'][st] = rr
    print(f"true varsigma={st:2d}: kappa bias {100*rr['kap_relbias']:6.2f}% sd {100*rr['kap_relstd']:.2f}% NMSE trk {rr['nmse_trk'][0]:.2f} orc {rr['nmse_orc'][0]:.2f} penalty {rr['penalty_dB']:.2f} dB")
rr = sweep([0.10], nr=300, model='saleh', model_tab='rapp')[0]; R['shape']['saleh_vs_rapp'] = rr
print(f"true Saleh, Rapp table: kappa abs bias {100*rr['kap_relbias']:.1f}% penalty {rr['penalty_dB']:.2f} dB")
json.dump(R, open('e_tracker.json', 'w'), default=float, indent=1)
