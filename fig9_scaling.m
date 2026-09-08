%FIG9_SCALING  Scaling in surface size N and compression smoothness (Fig. 9).
clear; rng(707); P = params(); P.nReal = 220;
PT = 10^(30/10); su2 = PT*P.rho_b + P.sig1_2;
[kap, sd2] = bussgang_coeff(su2, P.p, P.Asat, P.smooth, 'rapp', 2e6);
gam = sd2/(abs(kap)^2*su2);
fprintf('operating point: IBO=%.1f dB, gamma_d=%.2f dB\n\n', ...
        10*log10(P.Asat^2/(P.p^2*su2)), 10*log10(gam));

% ---- (a) sweep N with T = N+1 -------------------------------------------
Ns = [8 16 32 64]; RA = zeros(numel(Ns),3); TH = zeros(numel(Ns),1);
for i = 1:numel(Ns)
    o = struct('PT',PT,'Asat',P.Asat,'N',Ns(i));
    RA(i,1) = sim_nmse(P, setfield(setfield(o,'model','linear'),'gain',P.p));
    RA(i,2) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',P.p));
    RA(i,3) = sim_nmse(P, setfield(setfield(o,'model','rapp'),  'gain',kap));
    T = Ns(i)+1;
    TH(i) = (Ns(i)/T)*gam + P.sig2_2/(P.K*T*P.Asat^2*P.rho_g);   % eq. (10)
    fprintf('N=%3d  linear %7.2f conv %7.2f comp %7.2f | theory %7.2f\n', ...
            Ns(i), 10*log10(RA(i,:)), 10*log10(TH(i)));
end

% ---- (b) sweep the Rapp smoothness --------------------------------------
Sv = [0.5 1 2 4 8]; RB = zeros(numel(Sv),2); THB = zeros(numel(Sv),1);
for i = 1:numel(Sv)
    Q = P; Q.smooth = Sv(i);
    [k2, s2] = bussgang_coeff(su2, P.p, P.Asat, Sv(i), 'rapp', 1.5e6);
    o = struct('PT',PT,'Asat',P.Asat,'model','rapp');
    RB(i,1) = sim_nmse(Q, setfield(o,'gain',P.p));
    RB(i,2) = sim_nmse(Q, setfield(o,'gain',k2));
    g2 = s2/(abs(k2)^2*su2);
    THB(i) = (P.N/P.T)*g2 + P.sig2_2/(P.K*P.T*P.Asat^2*P.rho_g);
    fprintf('varsigma=%4.1f  |k|/p=%.4f  conv %7.2f comp %7.2f | theory %7.2f\n', ...
            Sv(i), abs(k2)/P.p, 10*log10(RB(i,:)), 10*log10(THB(i)));
end

figure;
subplot(1,2,1); semilogx(Ns,10*log10(RA(:,1)),'-o', Ns,10*log10(RA(:,3)),'-^', ...
    Ns,10*log10(TH),':', Ns,10*log10(RA(:,2)),'-s','LineWidth',1.8);
set(gca,'XTick',Ns); grid on; xlabel('RIS elements N (T=N+1)'); ylabel('NMSE [dB]');
legend('Linear','Compensated','Closed form','Conventional','Location','southwest');
title('(a) Scaling the surface');
subplot(1,2,2); semilogx(Sv,10*log10(RB(:,2)),'-^', Sv,10*log10(THB),':', ...
    Sv,10*log10(RB(:,1)),'-s','LineWidth',1.8);
set(gca,'XTick',Sv); grid on; xlabel('Rapp smoothness \varsigma'); ylabel('NMSE [dB]');
legend('Compensated','Closed form','Conventional','Location','southwest');
title('(b) Compression smoothness');
save('fig9_data.mat','Ns','RA','TH','Sv','RB','THB');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_scaling', 17);
