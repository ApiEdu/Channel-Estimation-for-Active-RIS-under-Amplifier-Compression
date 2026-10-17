"""Figs. 6-7 regenerated with the exact RMS-gain table + guard window (replaces the |kappa|/p version). seed 4242.
Fig 6: fixed calibration (nominal kappa) vs tracker vs oracle over A_sat deviation at P_T=30 dBm.
Fig 7: probe power ratio sweep at A_sat=0.14 (+3 dB)."""
import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import common as C
from common import P, db, cn
rng = np.random.default_rng(4242)
p, N, K = P['p'], P['N'], P['K']; PT2 = 10**3.0
xt, kt, rt = C.table('rapp', 2)
lg = np.log(xt)
def ratio_inv(rho):
    gi = lambda v: np.exp(np.interp(np.log(v), lg, np.log(rt)))
    ib = xt[xt*rho < xt[-1]]; rr = gi(ib)/gi(ib*rho); return ib, rr
def kappa_at(ib):
    return p*(np.interp(np.log(ib), lg, kt.real)+1j*np.interp(np.log(ib), lg, kt.imag))
def run(A_true, rho_dB, nr, guard=(-12, 8)):
    rho = 10**(rho_dB/10); ib_t, rr = ratio_inv(rho)
    su2 = PT2*P['rho_b']+P['sig1_2']; kap_true = C.moments('rapp', su2, p, A_true, 2)[0].real
    kap_nom = C.moments('rapp', su2, p, P['Asat'], 2)[0].real; ibo_true = C.ibo(A_true, p, su2)
    e = {'nom': [], 'trk': [], 'orc': []}; den = []; ibh = []
    for _ in range(nr):
        b = np.sqrt(P['rho_b'])*cn(rng, N); G = C.make_G(rng, K, N, P['rho_g']); hd = np.sqrt(P['rho_d'])*cn(rng, K)
        c2 = C.observe(b, G, hd, PT2, rng, 'rapp', A_true, 2); c1 = C.observe(b, G, hd, PT2/rho, rng, 'rapp', A_true, 2)
        r = np.linalg.norm(c2)/np.linalg.norm(c1)
        ib = float(np.exp(np.interp(r, rr, np.log(ib_t)))); ib = float(np.clip(ib, 10**(guard[0]/10), 10**(guard[1]/10)))
        kh = kappa_at(ib).real; ibh.append(ib)
        for k, kv in (('nom', kap_nom), ('trk', kh), ('orc', kap_true)):
            e[k].append(np.sum(np.abs(b-c2/kv)**2))
        den.append(np.sum(np.abs(b)**2))
    out = {k: C.boot_db(v, den, rng, 400) for k, v in e.items()}
    out['ibo_bias_dB'] = float(np.mean(db(ibh))-db(ibo_true)); out['ibo_true_dB'] = float(db(ibo_true))
    return out
R = {'dev': {}, 'probe': {}}
print('== Fig 6: A_sat deviation (P_T=30 dBm, varrho=10 dB, 300 realizations) ==')
for dev in (-6, -3, 0, 3, 6):
    A = P['Asat']*10**(dev/20); o = run(A, 10.0, 300); R['dev'][dev] = o
    print(f"dev {dev:+d} dB (IBO {o['ibo_true_dB']:.1f}): nominal {o['nom'][0]:6.2f} [{o['nom'][1]:.2f},{o['nom'][2]:.2f}] tracked {o['trk'][0]:6.2f} [{o['trk'][1]:.2f},{o['trk'][2]:.2f}] oracle {o['orc'][0]:6.2f}  tracked-oracle {o['trk'][0]-o['orc'][0]:+.2f}")
print('== Fig 7: probe ratio (A_sat=0.14) ==')
for v in (4, 7, 10, 14, 18, 25):
    o = run(0.14, v, 300); R['probe'][v] = o
    print(f"varrho {v:2d} dB: NMSE {o['trk'][0]:6.2f} [{o['trk'][1]:.2f},{o['trk'][2]:.2f}]  IBO bias {o['ibo_bias_dB']:+.2f} dB")
json.dump(R, open('e_tracker_main.json', 'w'), default=float, indent=1)
import plot_tracker_main  # draws Figs. 6 and 7 from e_tracker_main.json
