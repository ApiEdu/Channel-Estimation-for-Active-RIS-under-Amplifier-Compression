%FIG3_ANTENNA_SCALING  Antennas do not mitigate distortion (Fig. 3).
clear; rng(11); P = params(); P.nReal = 250;
PT = 10^(30/10); Ks = [1 2 4 8 16 32 64];
su2 = PT*P.rho_b + P.sig1_2;
kap = bussgang_coeff(su2, P.p, P.Asat, P.smooth, 'rapp', 2e6);
R = zeros(numel(Ks),3);
for i = 1:numel(Ks)
    o = struct('PT',PT,'Asat',P.Asat,'K',Ks(i));
    R(i,1) = sim_nmse(P, setfield(setfield(o,'model','linear'),'gain',P.p));
    R(i,2) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',P.p));
    R(i,3) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',kap));
    fprintf('K=%3d  %7.2f %7.2f %7.2f\n', Ks(i), 10*log10(R(i,:)));
end
fprintf('\nlinear gain K=1->64: %.1f dB | compensated: %.1f dB\n', ...
    10*log10(R(1,1))-10*log10(R(end,1)), 10*log10(R(1,3))-10*log10(R(end,3)));
figure; semilogx(Ks, 10*log10(R(:,1)),'-o', Ks, 10*log10(R(:,3)),'--^', ...
                 Ks, 10*log10(R(:,2)),'-s','LineWidth',1.9,'MarkerSize',6);
set(gca,'XTick',Ks,'XScale','log'); grid on;
xlabel('Number of BS antennas K'); ylabel('NMSE at P_T = 30 dBm  [dB]');
legend('Linear model','Compensated LS','Conventional LS','Location','southwest');
save('fig3_data.mat','Ks','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_antenna_scaling', 12);
