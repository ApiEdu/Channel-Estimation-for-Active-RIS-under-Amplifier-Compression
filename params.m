function P = params()
%PARAMS  Simulation parameters (Table 1 of the paper).
P.N      = 16;              % RIS elements
P.K      = 16;              % BS antennas
P.T      = P.N + 1;         % training slots
P.p      = 10;              % small-signal amplitude gain (20 dB)
P.Asat   = 0.10;            % saturation amplitude [sqrt(mW)]
P.smooth = 2;               % Rapp smoothness (varsigma)
P.rho_b  = 1e-6;            % user -> RIS large-scale gain
P.rho_g  = 1e-3/50^2;       % RIS  -> BS large-scale gain
P.rho_d  = 1e-7;            % direct link
P.sig1_2 = 10^(-70/10);     % RIS thermal noise [mW]
P.sig2_2 = 10^(-80/10);     % BS noise [mW]
P.nReal  = 600;             % Monte Carlo realizations
P.varrho = 10;              % probe power ratio [dB]
end
