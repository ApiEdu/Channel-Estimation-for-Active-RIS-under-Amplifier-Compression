"""Corrected error-floor theory (static bias + fluctuating part), K-scaling for scattering backward channels,
and Bayes-MMSE benchmark with the true Rapp model.  seed 7."""
import numpy as np, json
import common as C
from common import P, db, cn
rng = np.random.default_rng(7)
p, N, K = P['p'], P['N'], P['K']; A0 = P['Asat']

def decomp(PT, s1, M=4000, L=400, A=A0, s=2):
    """sigma_s^2 (static, b-dependent), sigma_f^2 (fluctuating incl. kappa*v), kappa, sigma_d^2 by nested Monte Carlo."""
    r = np.random.default_rng(99)
    b = np.sqrt(P['rho_b'])*cn(r, M); v = np.sqrt(s1)*cn(r, M, L)
    u = np.sqrt(PT)*b[:, None]+v; g = C.rapp(u, p, A, s)
    kap = np.mean(np.conj(u)*g)/(PT*P['rho_b']+s1)
    d = g-kap*u; dbar = d.mean(1)
    return dict(kap=float(kap.real), sd2=float(np.mean(np.abs(d)**2)), ss2=float(np.mean(np.abs(dbar)**2)),
                sf2=float(np.mean(np.abs(d-dbar[:, None])**2)+kap.real**2*s1))

def theory(PT, s1, T, Kk, geom='rank1', s2=P['sig2_2']):
    m = decomp(PT, s1); k2 = m['kap']**2; den = k2*PT*P['rho_b']
    if geom == 'rank1':
        fl = N/T*m['sf2']; rx = s2/(Kk*PT*T*k2*P['rho_g'])
    else:   # independent Rayleigh, MRC across antennas (K>=2)
        fl = (1+(N-1)/(Kk-1))*m['sf2']/T; rx = s2/(T*(Kk-1)*PT*k2*P['rho_g'])
    return (m['ss2']+fl)/den+rx/P['rho_b'], m

def sim_comp(PT, s1, T, Kk, geom='rank1', nr=300, Kr=None, mode='comp', kap=None):
    m = decomp(PT, s1); kap = m['kap']
    err = []; den = []
    for _ in range(nr):
        b = np.sqrt(P['rho_b'])*cn(rng, N); G = C.make_G(rng, Kk, N, P['rho_g'], geom, Kr); hd = np.sqrt(P['rho_d'])*cn(rng, Kk)
        ch = C.observe(b, G, hd, PT, rng, 'rapp', A0, 2, T=T, s1=s1)
        err.append(np.sum(np.abs(b-ch/kap)**2)); den.append(np.sum(np.abs(b)**2))
    return C.boot_db(err, den, rng, 400)

R = {}
print("== validation of corrected closed form (comp-LS NMSE, dB; sim [95% CI] | new | old (N/T)gamma_d) ==")
R['val'] = []
for (T, ptdb, s1db, Kk, geom) in [(17,30,-70,16,'rank1'),(34,30,-70,16,'rank1'),(68,30,-70,16,'rank1'),(136,30,-70,16,'rank1'),
                                   (17,30,-30,16,'rank1'),(34,30,-30,16,'rank1'),(68,30,-30,16,'rank1'),(136,30,-30,16,'rank1'),
                                   (17,30,-30,64,'rayleigh'),(68,30,-30,64,'rayleigh'),(17,30,-70,64,'rayleigh'),(17,20,-50,16,'rank1')]:
    PT = 10**(ptdb/10); s1 = 10**(s1db/10)
    sm = sim_comp(PT, s1, T, Kk, geom)
    th, m = theory(PT, s1, T, Kk, geom)
    old = (N/T*m['sd2']/(m['kap']**2*PT*P['rho_b']) + P['sig2_2']/(Kk*PT*T*m['kap']**2*P['rho_g'])/P['rho_b'])
    row = dict(T=T, PT_dBm=ptdb, s1_dBm=s1db, K=Kk, geom=geom, sim=sm, new=float(db(th)), old=float(db(old)),
               gs=float(db(m['ss2']/(m['kap']**2*PT*P['rho_b']))), gf=float(db(m['sf2']/(m['kap']**2*PT*P['rho_b']))))
    R['val'].append(row)
    print(f"T={T:3d} P={ptdb}dBm s1={s1db}dBm K={Kk:3d} {geom:8s} sim {sm[0]:6.2f} [{sm[1]:.2f},{sm[2]:.2f}] | new {row['new']:6.2f} | old {row['old']:6.2f}  (gs {row['gs']:.1f}, gf {row['gf']:.1f} dB)")

print("\n== K-scaling for four backward-channel geometries (comp-LS, P_T=30 dBm, s1=-70 dBm) ==")
PT = 1e3; s1 = P['sig1_2']; T = N+1
R['K'] = {}
for lab, geom, Kr in (('rank-one LoS','rank1',None), ('Rician 10 dB','rician',10), ('Rician 0 dB','rician',0), ('Rayleigh','rayleigh',None)):
    rows = []
    for Kk in (2, 4, 8, 16, 32, 64, 128):
        sm = sim_comp(PT, s1, T, Kk, geom, nr=250, Kr=Kr)
        th = theory(PT, s1, T, Kk, 'rank1' if geom == 'rank1' else 'rayleigh')[0] if geom in ('rank1', 'rayleigh') else None
        rows.append(dict(K=Kk, sim=sm, th=None if th is None else float(db(th))))
    R['K'][lab] = rows
    print(lab.ljust(13), ' '.join(f"K={r['K']}:{r['sim'][0]:6.2f}" for r in rows), '| theory:', ' '.join('--' if r['th'] is None else f"{r['th']:.2f}" for r in rows))
# linear reference (no compression): NMSE ~ rx term only
lin = []
for Kk in (2, 4, 8, 16, 32, 64, 128):
    err = []; den = []
    for _ in range(250):
        b = np.sqrt(P['rho_b'])*cn(rng, N); G = C.make_G(rng, Kk, N, P['rho_g'], 'rank1'); hd = np.sqrt(P['rho_d'])*cn(rng, Kk)
        ch = C.observe(b, G, hd, PT, rng, 'linear', None, 2, T=T)
        err.append(np.sum(np.abs(b-ch/p)**2)); den.append(np.sum(np.abs(b)**2))
    lin.append(C.boot_db(err, den, rng, 400))
R['K']['linear'] = lin
print('linear model  ', ' '.join(f"K={k}:{l[0]:6.2f}" for k, l in zip((2,4,8,16,32,64,128), lin)))

print("\n== Bayes-MMSE benchmark with the true Rapp model (genie noise level), T=N+1, K=16 ==")
R['bayes'] = []
bs_prior = np.sqrt(P['rho_b'])*cn(rng, 40000)
for ptdb in (0, 5, 10, 15, 20, 25, 30, 35, 40, 45):
    PT = 10**(ptdb/10); s1 = P['sig1_2']; m = decomp(PT, s1); kap = m['kap']
    gs = C.rapp(np.sqrt(PT)*bs_prior, p, A0, 2)/np.sqrt(PT)
    CH = []; B = []
    for _ in range(300):
        b = np.sqrt(P['rho_b'])*cn(rng, N); G = C.make_G(rng, K, N, P['rho_g']); hd = np.sqrt(P['rho_d'])*cn(rng, K)
        CH.append(C.observe(b, G, hd, PT, rng, 'rapp', A0, 2)); B.append(b)
    CH = np.array(CH); B = np.array(B)
    nu = np.mean(np.abs(CH-C.rapp(np.sqrt(PT)*B, p, A0, 2)/np.sqrt(PT))**2)
    den = np.sum(np.abs(B)**2, 1)
    e_conv = np.sum(np.abs(B-CH/p)**2, 1); e_comp = np.sum(np.abs(B-CH/kap)**2, 1)
    alpha = np.sum(np.real(np.conj(B)*CH/kap))/np.sum(np.abs(CH/kap)**2)
    e_bl = np.sum(np.abs(B-alpha*CH/kap)**2, 1)
    e_mm = np.zeros(300)
    for i in range(300):
        for n in range(N):
            w = np.exp(-np.abs(CH[i, n]-gs)**2/nu); w /= w.sum()
            e_mm[i] += np.abs(B[i, n]-(w*bs_prior).sum())**2
    row = dict(PT_dBm=ptdb, conv=C.boot_db(e_conv, den, rng, 400), comp=C.boot_db(e_comp, den, rng, 400),
               blmmse=C.boot_db(e_bl, den, rng, 400), bayes=C.boot_db(e_mm, den, rng, 400))
    R['bayes'].append(row)
    print(f"P={ptdb:2d} dBm conv {row['conv'][0]:6.2f} comp {row['comp'][0]:6.2f} BLMMSE {row['blmmse'][0]:6.2f} [{row['blmmse'][1]:.2f},{row['blmmse'][2]:.2f}] Bayes-MMSE {row['bayes'][0]:6.2f} [{row['bayes'][1]:.2f},{row['bayes'][2]:.2f}]  gap BL-Bayes {row['blmmse'][0]-row['bayes'][0]:.2f} dB")
json.dump(R, open('e_theory.json', 'w'), default=float, indent=1)
