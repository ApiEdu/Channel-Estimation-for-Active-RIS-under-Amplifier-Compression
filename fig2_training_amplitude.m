%FIG2_TRAINING_AMPLITUDE  The prior-art training rule (Fig. 2).
%   Sum-NMSE over [h_d ; b] is used because the amplitude trade-off in the
%   prior art originates in the direct link.
clear; rng(31); P = params(); P.nReal = 700;
PT = 10^(30/10);
A  = [1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8 9 10];
lin = zeros(size(A)); nl = zeros(size(A));
for i = 1:numel(A)
    o = struct('PT',PT,'Asat',P.Asat,'a',A(i),'gain',A(i),'sumMSE',true);
    lin(i) = sim_nmse(P, setfield(o,'model','linear'));
    nl(i)  = sim_nmse(P, setfield(o,'model','rapp'));
    fprintf('a=%5.2f  linear %7.2f  compressed %7.2f\n', ...
            A(i), 10*log10(lin(i)), 10*log10(nl(i)));
end
[~,im] = min(nl);
fprintf('\ntrue optimum a=%.1f; penalty at a=%.0f is %.2f dB\n', ...
        A(im), A(end), 10*log10(nl(end))-10*log10(nl(im)));
figure; hold on; grid on; box on;
plot(A, 10*log10(lin), '-o','LineWidth',1.9,'MarkerSize',6);
plot(A, 10*log10(nl),  '-s','LineWidth',1.9,'MarkerSize',6);
xlabel('Training amplification a'); ylabel('Sum-NMSE  [dB]');
legend('Linear model (prior-art assumption)','With amplifier compression', ...
       'Location','southwest');
save('fig2_data.mat','A','lin','nl');

% --- export for the manuscript (MDPI: >=300 dpi) ---
mdpi_style('fig_training_amplitude', 12);
