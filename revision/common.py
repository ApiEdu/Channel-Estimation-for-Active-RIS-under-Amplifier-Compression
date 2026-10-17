"""Shared code for the revision experiments (Telecom revision, Oct 2026).
Mirrors matlab/params.m, rapp.m, saleh.m, sim_nmse.m.  Python 3 + NumPy/SciPy/Matplotlib.
Every script fixes its RNG seed; Monte Carlo CIs are 95 % bootstrap intervals over realizations."""
import numpy as np
db = lambda x: 10*np.log10(x)
trap = getattr(np, 'trapezoid', None) or np.trapz

P = dict(N=16, K=16, p=10.0, Asat=0.10, smooth=2, rho_b=1e-6, rho_g=1e-3/50**2,
         rho_d=1e-7, sig1_2=10**(-70/10), sig2_2=10**(-80/10), varrho_dB=10.0)

# ---------------------------------------------------------------- amplifiers
def rapp(u, p, A, s=2):
    return p*u/(1+(p*np.abs(u)/A)**(2*s))**(1/(2*s))

def saleh(u, p, A, bA=1.1517, aP=4.0033, bP=9.1040):
    us = 2*np.sqrt(bA)*A/p
    x = np.abs(u)/us
    return p*np.abs(u)/(1+bA*x**2)*np.exp(1j*(np.angle(u)+aP*x**2/(1+bP*x**2)))

def amp(model, u, p, A, s=2):
    return rapp(u, p, A, s) if model == 'rapp' else saleh(u, p, A)

# -------------------------------------------------- Bussgang moments (quadrature)
def moments(model, sig2, p, A, s=2, n=60001):
    """kappa (complex), sigma_d^2, rms gain^2 = E|g|^2/sig2 for u ~ CN(0,sig2).
    t=|u|^2/sig2 ~ Exp(1);  E[...] = int e^-t (...) dt  (phase of u cancels)."""
    t = np.linspace(0, 40, n); w = np.exp(-t)
    r = np.sqrt(t*sig2)
    g = amp(model, r.astype(complex), p, A, s)           # u = r (real positive)
    G = np.where(r > 0, g/np.where(r > 0, r, 1), p)       # complex gain g/u
    kap = trap(w*t*G, t)
    pw = trap(w*t*np.abs(G)**2, t)                        # E|g|^2/sig2
    return kap, (pw-abs(kap)**2)*sig2, pw

def ibo(A, p, sig2): return A**2/(p**2*sig2)

def table(model='rapp', s=2, p=10.0, lo=-5, hi=5, npts=201):
    """Universal curves versus IBO (sigma_u^2 normalised to 1, A = sqrt(IBO) p)."""
    x = np.logspace(lo, hi, npts)
    kap = np.empty(npts, complex); rms = np.empty(npts)
    for i, v in enumerate(x):
        k, _, pw = moments(model, 1.0, p, np.sqrt(v)*p, s)
        kap[i] = k/p; rms[i] = np.sqrt(pw)/p
    return x, kap, rms

# ------------------------------------------------------------------ simulator
def cn(rng, *sh): return (rng.standard_normal(sh)+1j*rng.standard_normal(sh))/np.sqrt(2)

def dft(T, N): return np.exp(-2j*np.pi*np.outer(np.arange(T), np.arange(1, N+1))/T)

def make_G(rng, K, N, rho_g, geom='rank1', Kr_dB=None):
    if geom == 'rank1':
        return np.sqrt(rho_g)*np.outer(np.exp(1j*np.pi*np.arange(K)*np.sin(rng.uniform(-1.5, 1.5))),
                                        np.conj(np.exp(1j*np.pi*np.arange(N)*np.sin(rng.uniform(-1.5, 1.5)))))
    if geom == 'rayleigh':
        return np.sqrt(rho_g)*cn(rng, K, N)
    if geom == 'rician':
        k = 10**(Kr_dB/10); a = np.sqrt(k/(k+1))
        los = make_G(rng, K, N, 1.0, 'rank1')
        return np.sqrt(rho_g)*(a*los+np.sqrt(1-a*a)*cn(rng, K, N))
    raise ValueError(geom)

def observe(b, G, hd, PT, rng, model='rapp', A=None, s=2, T=None, s1=None, s2=None, p=None, a=None):
    """One training run.  Returns chat = MRC estimate of E_t[g(u)]/sqrt(PT)  (so b_conv = chat/p, b_comp = chat/kappa).
    A may be scalar or length-N vector; a = training gain (defaults to p)."""
    N = b.size; K = G.shape[0]
    T = T or N+1; s1 = P['sig1_2'] if s1 is None else s1; s2 = P['sig2_2'] if s2 is None else s2
    p = P['p'] if p is None else p; A = P['Asat'] if A is None else A
    F = dft(T, N)
    v = np.sqrt(s1)*cn(rng, T, N); z = np.sqrt(s2)*cn(rng, T, K)
    u = np.sqrt(PT)*b[None, :]+v
    out = amp(model, u, p if a is None else a, A, s) if model != 'linear' else (p if a is None else a)*u
    Y = (out*F)@G.T+np.sqrt(PT)*hd[None, :]+z
    proj = F.conj().T@Y/T
    c = (np.conj(G).T*proj).sum(1)/(np.abs(G)**2).sum(0)
    return c/np.sqrt(PT)

def boot_db(err, den, rng, B=1000):
    """NMSE = sum(err)/sum(den) with 95 % bootstrap CI (dB)."""
    err = np.asarray(err); den = np.asarray(den); n = err.size
    idx = rng.integers(0, n, (B, n))
    bs = db(err[idx].sum(1)/den[idx].sum(1))
    return db(err.sum()/den.sum()), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
