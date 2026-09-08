%FIG1_NMSE_FLOOR  Error floor and non-monotonicity (Fig. 1).
clear; rng(7); P = params();
PT_dBm = -10:2.5:45;  PT = 10.^(PT_dBm/10);
lab = {'linear','conv','comp','blmmse'};
R = zeros(numel(PT), 4);  IBO = zeros(size(PT));

for i = 1:numel(PT)
    su2 = PT(i)*P.rho_b + P.sig1_2;
    [kap, sd2] = bussgang_coeff(su2, P.p, P.Asat, P.smooth, 'rapp', 1e6);
    IBO(i) = P.Asat^2/(P.p^2*su2);
    ak2 = abs(kap)^2;
    eps_b = P.N*(ak2*P.sig1_2 + sd2)/(PT(i)*P.T*ak2) ...
          + P.sig2_2/(P.K*PT(i)*P.T*ak2*P.rho_g);         % eq. (9)
    shrink = P.rho_b/(P.rho_b + eps_b);                    % BLMMSE, eq. (12)
    o = struct('PT',PT(i),'Asat',P.Asat);
    R(i,1) = sim_nmse(P, setfield(setfield(o,'model','linear'),'gain',P.p));
    R(i,2) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',P.p));
    R(i,3) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',kap));
    oo = setfield(setfield(o,'model','rapp'),'gain',kap); oo.shrink = shrink;
    R(i,4) = sim_nmse(P, oo);
    fprintf('PT=%5.1f dBm  IBO=%6.1f dB  %7.2f %7.2f %7.2f %7.2f\n', ...
        PT_dBm(i), 10*log10(IBO(i)), 10*log10(R(i,:)));
end

figure; hold on; grid on; box on;
plot(PT_dBm, 10*log10(R(:,1)), '-',   'LineWidth',1.9);
plot(PT_dBm, 10*log10(R(:,2)), '-s',  'LineWidth',1.7,'MarkerSize',5);
plot(PT_dBm, 10*log10(R(:,3)), '--^', 'LineWidth',1.7,'MarkerSize',5);
plot(PT_dBm, 10*log10(R(:,4)), '-.o', 'LineWidth',1.9,'MarkerSize',5);
xlabel('Pilot transmit power P_T  [dBm]'); ylabel('NMSE  [dB]');
legend('Linear model','Conventional LS','Compensated LS','Bussgang-LMMSE', ...
       'Location','southwest');
ylim([-46 14]);
save('fig1_data.mat','PT_dBm','IBO','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_nmse_floor_ibo', 12);
