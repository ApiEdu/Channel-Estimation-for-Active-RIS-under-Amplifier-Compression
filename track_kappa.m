function [kappa_hat, ibo_hat] = track_kappa(P, PT2, varrho_dB, Asat_true, model)
%TRACK_KAPPA  Algorithm 1: dual-power probing.
%   Estimates the operating input back-off from the ratio of two channel
%   estimates taken at pilot powers PT2 and PT2/varrho, then returns the
%   corresponding Bussgang gain. Requires neither Asat nor rho_b.
if nargin < 5, model = 'rapp'; end
varrho = 10^(varrho_dB/10);

% ---- offline: universal compression curve f(IBO) = |kappa|/p -------------
persistent ibo_grid f_grid cached
if isempty(cached) || ~strcmp(cached, model)
    ibo_grid = logspace(-3, 3, 60);
    f_grid   = zeros(size(ibo_grid));
    for i = 1:numel(ibo_grid)
        Aeq = sqrt(ibo_grid(i) * P.p^2 * 1.0);       % sigma_u2 normalised to 1
        k   = bussgang_coeff(1.0, P.p, Aeq, P.smooth, model, 2e5);
        f_grid(i) = abs(k)/P.p;
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
Lopt   = fminbnd(cost, log(1e-3), log(1e3));
ibo_hat   = exp(Lopt);
kappa_hat = P.p * f(ibo_hat);
end
