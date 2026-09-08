%FIG6_PROBE_DESIGN  Probe power ratio: observability vs probe noise (Fig. 6).
clear; rng(202); P = params(); P.nReal = 350;
PT = 10^(30/10); Asat_true = 0.14;          % +3 dB, where fixed calibration fails
V = [4 7 10 14 18 25];
su2 = PT*P.rho_b + P.sig1_2;
ibo_true = Asat_true^2/(P.p^2*su2);
nm = zeros(size(V)); bias = zeros(size(V));
for i = 1:numel(V)
    [kh, ih] = track_kappa(P, PT, V(i), Asat_true, 'rapp');
    o = struct('PT',PT,'model','rapp','Asat',Asat_true,'gain',kh);
    nm(i)   = sim_nmse(P, o);
    bias(i) = 10*log10(ih) - 10*log10(ibo_true);
    fprintf('varrho=%2d dB  NMSE %7.2f  IBO bias %+6.2f dB\n', ...
            V(i), 10*log10(nm(i)), bias(i));
end
[~,ib] = min(nm);
fprintf('\noptimum probe ratio: %d dB\n', V(ib));
figure;
yyaxis left;  plot(V, 10*log10(nm), '-o','LineWidth',2.0,'MarkerSize',7);
ylabel('NMSE  [dB]');
yyaxis right; plot(V, bias, '--s','LineWidth',1.6,'MarkerSize',5);
ylabel('Bias of estimated IBO  [dB]');
xlabel('Probe power ratio \varrho = P_2/P_1  [dB]'); grid on;
save('fig6_data.mat','V','nm','bias');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_probe_design', 12);
