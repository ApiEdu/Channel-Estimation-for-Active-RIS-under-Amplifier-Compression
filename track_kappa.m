function [kappa_hat, ibo_hat] = track_kappa(P, PT2, varrho_dB, Asat_true, model)
%TRACK_KAPPA  Algorithm 1: dual-power probing.
%   Estimates the operating input back-off from the ratio of two channel
%   estimates taken at pilot powers PT2 and PT2/varrho, then returns the
%   corresponding Bussgang gain. Requires neither Asat nor rho_b.
if nargin < 5, model = 'rapp'; end
varrho = 10^(varrho_dB/10);

% ---- offline: universal compression curve h(IBO) = RMS gain ---------------
persistent ibo_grid f_grid k_grid cached
if isempty(cached) || ~strcmp(cached, model)
    ibo_grid = logspace(-3, 3, 60);
    f_grid   = zeros(size(ibo_grid));      % RMS gain h(IBO), eq. (14)
    k_grid   = zeros(size(ibo_grid));      % complex kappa(IBO)
    for i = 1:numel(ibo_grid)
        Aeq = sqrt(ibo_grid(i) * P.p^2 * 1.0);       % sigma_u2 normalised to 1
        [k, sd2] = bussgang_coeff(1.0, P.p, Aeq, P.smooth, model, 2e5);
        f_grid(i) = sqrt(abs(k)^2 + sd2)/P.p;     % RMS gain, not |kappa|/p
        k_grid(i) = k;
    end
    cached = model;
end
f = @(x) interp1(log(ibo_grid), f_grid, log(min(max(x,1e-3),1e3)), 'pchip');

% ---- probe ---------------------------------------------------------------
o = struct('model',model,'Asat',Asat_true,'gain',P.p,'a',P.p);
o.PT = PT2;         [~, n2] = sim_nmse(P, o);
o.PT = PT2/varrho;  [~, n1] = sim_nmse(P, o);
r = n2/n1;

% ---- solve f(IBO)/f(varrho*IBO) = r --------------------------------------
cost   = @(L) (f(exp(L))/f(varrho*exp(L)) - r)^2;
% guard window [-12,+8] dB: outside it the ratio is flat and the root ill-conditioned
Lopt   = fminbnd(cost, log(10^(-12/10)), log(10^(8/10)));
ibo_hat   = exp(Lopt);
% k_grid already holds kappa itself (it includes the gain p): do NOT multiply by p again.
% Real and imaginary parts are interpolated separately (kappa is complex for Saleh).
kappa_hat = interp1(log(ibo_grid), real(k_grid), log(ibo_hat), 'pchip') ...
          + 1j*interp1(log(ibo_grid), imag(k_grid), log(ibo_hat), 'pchip');
end
