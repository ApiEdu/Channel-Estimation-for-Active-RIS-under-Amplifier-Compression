%FIG7_SALEH  AM/PM conversion and the validity limit of compensation (Fig. 7).
clear; rng(404); P = params(); P.nReal = 300;
PT_dBm = [10 20 25 30 35 40];
R = zeros(numel(PT_dBm),3); IBO = zeros(numel(PT_dBm),1); phi = IBO;
for i = 1:numel(PT_dBm)
    PT  = 10^(PT_dBm(i)/10);
    su2 = PT*P.rho_b + P.sig1_2;
    kap = bussgang_coeff(su2, P.p, P.Asat, [], 'saleh', 1.5e6);   % COMPLEX
    IBO(i) = 10*log10(P.Asat^2/(P.p^2*su2));
    phi(i) = rad2deg(angle(kap));
    o = struct('PT',PT,'model','saleh','Asat',P.Asat);
    R(i,1) = sim_nmse(P, setfield(o,'gain',P.p));        % conventional
    R(i,2) = sim_nmse(P, setfield(o,'gain',abs(kap)));   % magnitude only
    R(i,3) = sim_nmse(P, setfield(o,'gain',kap));        % complex kappa
    fprintf('IBO=%6.1f dB  arg(k)=%5.1f deg | conv %7.2f  |k| %7.2f  complex %7.2f\n', ...
            IBO(i), phi(i), 10*log10(R(i,:)));
end
figure; hold on; grid on; box on;
plot(IBO, 10*log10(R(:,1)),'-s','LineWidth',1.8,'MarkerSize',6);
plot(IBO, 10*log10(R(:,2)),'-.d','LineWidth',1.8,'MarkerSize',6);
plot(IBO, 10*log10(R(:,3)),'-o','LineWidth',2.1,'MarkerSize',7);
set(gca,'XDir','reverse');
xlabel('Input back-off  [dB]'); ylabel('NMSE  [dB]');
legend('Conventional LS','Compensated, |\kappa| only','Compensated, complex \kappa', ...
       'Location','southwest');
save('fig7_data.mat','IBO','phi','R');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_saleh', 12);
