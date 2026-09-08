%FIG4_GAIN_MISMATCH  Calibration asymmetry (Fig. 4).
%   x-axis is (true - assumed) in dB. POSITIVE means the true Asat is LARGER
%   than assumed, i.e. we assumed MORE compression than occurs -> kappa too
%   small -> the estimate is shrunk -> benign. NEGATIVE is the dangerous side.
clear; rng(11); P = params(); P.nReal = 800;
PT = 10^(30/10); Asat_nom = P.Asat;
mm = -10:1:10;  Asat_true = Asat_nom*10.^(mm/20);
su2 = PT*P.rho_b + P.sig1_2;
kap_nom = bussgang_coeff(su2, P.p, Asat_nom, P.smooth, 'rapp', 2e6);
R = zeros(numel(mm),3);
for i = 1:numel(mm)
    kap_true = bussgang_coeff(su2, P.p, Asat_true(i), P.smooth, 'rapp', 2e6);
    o = struct('PT',PT,'model','rapp','Asat',Asat_true(i));
    R(i,1) = sim_nmse(P, setfield(o,'gain',P.p));       % no compensation
    R(i,2) = sim_nmse(P, setfield(o,'gain',kap_nom));   % fixed calibration
    R(i,3) = sim_nmse(P, setfield(o,'gain',kap_true));  % oracle
    fprintf('mismatch %+5.1f dB  naive %7.2f  fixed %7.2f  oracle %7.2f\n', ...
            mm(i), 10*log10(R(i,:)));
end
figure; hold on; grid on; box on;
plot(mm, 10*log10(R(:,1)),'-s','LineWidth',1.7,'MarkerSize',5);
plot(mm, 10*log10(R(:,2)),'-o','LineWidth',1.9,'MarkerSize',6);
plot(mm, 10*log10(R(:,3)),'--^','LineWidth',1.7,'MarkerSize',5);
xlabel('A_{sat} mismatch: true - assumed  [dB]');
ylabel('NMSE at P_T = 30 dBm  [dB]');
legend('Conventional LS','Compensated, nominal calibration','Oracle', ...
       'Location','northwest');
save('fig4_data.mat','mm','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_gain_mismatch', 12);
