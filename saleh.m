function y = saleh(u, p, Asat)
% SALEH  AM/AM and AM/PM amplifier, eq. (13).
%   u_s is set so the PEAK output amplitude equals Asat, matching the Rapp
%   convention. Because of the AM/PM term, kappa is COMPLEX.
bA = 1.1517; aP = 4.0033; bP = 9.1040;   % Saleh TWTA parameters
u_s = 2*sqrt(bA)*Asat/p;
x   = abs(u)/u_s;
amp = p*abs(u) ./ (1 + bA*x.^2);
ph  = aP*x.^2 ./ (1 + bP*x.^2);
y   = amp .* exp(1j*(angle(u) + ph));
end
