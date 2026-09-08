%FIG5_KAPPA_TRACKING  Online tracking vs fixed calibration (Fig. 5).
clear; rng(101); P = params(); P.nReal = 120;
PT = 10^(30/10); Asat_set = [0.05 0.07 0.10 0.14 0.20];
su2 = PT*P.rho_b + P.sig1_2;
kap_nom = bussgang_coeff(su2, P.p, 0.10, P.smooth, 'rapp', 2e6);
R = zeros(numel(Asat_set),3);
for i = 1:numel(Asat_set)
    At = Asat_set(i);
    kap_true = bussgang_coeff(su2, P.p, At, P.smooth, 'rapp', 2e6);
    kap_hat  = track_kappa(P, PT, P.varrho, At, 'rapp');   % Algorithm 1
    o = struct('PT',PT,'model','rapp','Asat',At);
    R(i,1) = sim_nmse(P, setfield(o,'gain',kap_nom));
    R(i,2) = sim_nmse(P, setfield(o,'gain',kap_hat));
    R(i,3) = sim_nmse(P, setfield(o,'gain',kap_true));
    fprintf('Asat=%.2f  k_true=%.3f k_hat=%.3f | nominal %7.2f tracked %7.2f oracle %7.2f\n', ...
            At, abs(kap_true), abs(kap_hat), 10*log10(R(i,:)));
end
dev = 20*log10(Asat_set/0.10);
figure; hold on; grid on; box on;
plot(dev, 10*log10(R(:,1)),'-s','LineWidth',1.8,'MarkerSize',6);
plot(dev, 10*log10(R(:,2)),'-o','LineWidth',2.1,'MarkerSize',7);
plot(dev, 10*log10(R(:,3)),'--^','LineWidth',1.7,'MarkerSize',6);
xlabel('A_{sat} deviation from nominal  [dB]');
ylabel('NMSE at P_T = 30 dBm  [dB]');
legend('Fixed calibration','Proposed dual-power tracking','Oracle \kappa', ...
       'Location','southwest');
save('fig5_data.mat','dev','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_kappa_tracking', 12);
