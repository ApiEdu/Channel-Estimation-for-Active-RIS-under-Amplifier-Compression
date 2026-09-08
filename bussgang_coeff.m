function [kappa, sigma_d2] = bussgang_coeff(sigma_u2, p, Asat, s, model, nMC)
% BUSSGANG_COEFF  Numerical Bussgang gain and distortion power, eq. (4).
%   kappa    = E{u* g(u)} / E{|u|^2}          (complex for the Saleh model)
%   sigma_d2 = E{|g(u)|^2} - |kappa|^2 sigma_u2
%   model : 'rapp' (default) or 'saleh'
if nargin < 5 || isempty(model), model = 'rapp'; end
if nargin < 6 || isempty(nMC),   nMC   = 4e5;    end
u = sqrt(sigma_u2/2)*(randn(nMC,1) + 1j*randn(nMC,1));
if strcmpi(model,'saleh'), g = saleh(u, p, Asat);
else,                      g = rapp(u, p, Asat, s); end
Pu       = mean(abs(u).^2);
kappa    = mean(conj(u).*g) / Pu;
sigma_d2 = max(mean(abs(g).^2) - abs(kappa)^2*Pu, 0);
end
