%FIG8_ELEMENT_SPREAD  Per-element saturation spread (Fig. 8).
%   The per-element oracle is flat: the spread costs nothing physically. The
%   loss comes only from approximating {kappa_n} by a single scalar.
clear; rng(606); P = params(); P.nReal = 300;
PT = 10^(30/10); su2 = PT*P.rho_b + P.sig1_2;
sigA = [0 0.5 1 2 3 4];

% tabulate kappa(Asat) once
grid_A = logspace(log10(0.02), log10(0.6), 40);
kg = zeros(size(grid_A));
for i = 1:numel(grid_A)
    kg(i) = bussgang_coeff(su2, P.p, grid_A(i), P.smooth, 'rapp', 3e5);
end
kap_of = @(a) interp1(grid_A, kg, a, 'pchip');

R = zeros(numel(sigA),4);
for i = 1:numel(sigA)
    As = P.Asat * 10.^(randn(1,P.N)*sigA(i)/20);       % per-element spread
    o = struct('PT',PT,'model','rapp','Asat',As);
    R(i,1) = sim_nmse(P, setfield(o,'gain',P.p));                 % conventional
    R(i,2) = sim_nmse(P, setfield(o,'gain',kap_of(P.Asat)));      % nominal kappa
    R(i,3) = sim_nmse(P, setfield(o,'gain',mean(kap_of(As))));    % mean kappa
    R(i,4) = sim_nmse(P, setfield(o,'gain',kap_of(As)));          % per element
    fprintf('sigma_A=%3.1f dB | conv %7.2f nominal %7.2f mean %7.2f per-element %7.2f\n', ...
            sigA(i), 10*log10(R(i,:)));
end
figure; hold on; grid on; box on;
plot(sigA, 10*log10(R(:,4)),'--^','LineWidth',1.8,'MarkerSize',6);
plot(sigA, 10*log10(R(:,3)),'-o', 'LineWidth',2.0,'MarkerSize',6);
plot(sigA, 10*log10(R(:,2)),'-d', 'LineWidth',1.8,'MarkerSize',6);
plot(sigA, 10*log10(R(:,1)),'-s', 'LineWidth',1.8,'MarkerSize',6);
xlabel('Per-element A_{sat} spread \sigma_A  [dB]');
ylabel('NMSE at P_T = 30 dBm  [dB]');
legend('Per-element \kappa_n (oracle)','Single \kappa, mean', ...
       'Single \kappa, nominal','Conventional LS','Location','southwest');
save('fig8_data.mat','sigA','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_element_spread', 12);
