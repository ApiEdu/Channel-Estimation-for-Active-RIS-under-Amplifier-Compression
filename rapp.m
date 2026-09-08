function y = rapp(u, p, Asat, s)
% RAPP  Soft-limiting reflection amplifier, eq. (2).
%   As Asat -> Inf this reduces exactly to y = p*u, i.e. the conventional
%   linear active-RIS model. AM/AM only: kappa comes out real.
if nargin < 4, s = 2; end
y = (p*u) ./ (1 + (p*abs(u)/Asat).^(2*s)).^(1/(2*s));
end
