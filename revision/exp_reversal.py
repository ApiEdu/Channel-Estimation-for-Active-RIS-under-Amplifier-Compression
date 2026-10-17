"""When does the maximum-gain training rule reverse?  a* (argmin sum-NMSE) versus direct-link strength,
amplifier headroom and receiver noise.  Mirrors matlab/fig2 (sum-NMSE over [h_d;b], estimator assumes gain a).  seed 31."""
import numpy as np, json
import common as C
from common import P, db, cn
rng = np.random.default_rng(31)
N, K, p = P['N'], P['K'], P['p']; T = N+1; F = C.dft(T, N); PT = 10**3.0
A_GRID = np.array([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10.0])

def sumnmse(a, model, rho_d, A, s2, nr=120):
    num = den = 0.0
    for _ in range(nr):
        b = np.sqrt(P['rho_b'])*cn(rng, N); hd = np.sqrt(rho_d)*cn(rng, K)
        G = C.make_G(rng, K, N, P['rho_g'])
        v = np.sqrt(P['sig1_2'])*cn(rng, T, N); z = np.sqrt(s2)*cn(rng, T, K)
        u = np.sqrt(PT)*b[None, :]+v
        out = a*u if model == 'linear' else C.rapp(u, a, A, 2)
        Y = (out*F)@G.T+np.sqrt(PT)*hd[None, :]+z
        h_est = Y.mean(0)/np.sqrt(PT)                           # DC bin (F columns orthogonal to DC)
        proj = F.conj().T@Y/T                                    # (N,K)
        b_est = proj/(a*G.T*np.sqrt(PT))                         # per-antenna LS, assumed gain a
        num += np.sum(np.abs(hd-h_est)**2)+np.sum(np.abs(b[:, None]-b_est)**2)/K
        den += np.sum(np.abs(hd)**2)+np.sum(np.abs(b)**2)
    return num/den

def best(model, **kw):
    v = np.array([sumnmse(a, model, **kw) for a in A_GRID])
    i = int(np.argmin(v)); return A_GRID[i], db(v[i]), db(v[-1])-db(v[i])

R = {}
print("== baseline (paper settings) ==")
base = dict(rho_d=P['rho_d'], A=P['Asat'], s2=P['sig2_2'])
b_nl = best('rapp', **base); b_lin = best('linear', **base)
print("compressed: a*=%.1f  penalty of a=10: %.2f dB | linear: a*=%.1f penalty %.2f dB" % (b_nl[0], b_nl[2], b_lin[0], b_lin[2]))
R['base'] = dict(nl=b_nl, lin=b_lin)
def scan(name, key, vals, fmt):
    R[name] = []
    print("\n== %s ==" % name)
    for v in vals:
        kw = dict(base); kw[key] = v
        nl = best('rapp', **kw); li = best('linear', **kw)
        R[name].append(dict(x=float(v), a_star=float(nl[0]), pen_amax=float(nl[2]), a_star_lin=float(li[0]), pen_lin=float(li[2])))
        print(fmt(v) + " | compressed a*=%4.1f penalty(a=10)=%6.2f dB | linear a*=%4.1f penalty %5.2f dB" % (nl[0], nl[2], li[0], li[2]))
scan('direct_link', 'rho_d', 10**np.arange(-9.0, -4.9, 1.0), lambda v: f"rho_d/rho_b = {10*np.log10(v/P['rho_b']):6.1f} dB")
scan('headroom', 'A', P['Asat']*10**(np.array([-6, -3, 0, 3, 6, 10, 14])/20), lambda v: f"Asat = {20*np.log10(v/P['Asat']):+5.0f} dB vs nominal")
scan('rx_noise', 's2', 10**(np.array([-110, -100, -90, -80, -70, -60])/10), lambda v: f"sigma2^2 = {10*np.log10(v):6.0f} dBm")
json.dump(R, open('e_reversal.json', 'w'), default=float, indent=1)
