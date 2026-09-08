function [nmse, bhat_norm] = sim_nmse(P, opt)
%SIM_NMSE  Monte Carlo NMSE of the user--RIS channel estimate.
%
%   opt fields (all optional, defaults in brackets):
%     .model    'rapp' | 'saleh' | 'linear'      ['rapp']
%     .PT       pilot power [mW]                 [1000]
%     .gain     gain assumed by the estimator    [P.p]   (complex allowed)
%     .Asat     saturation, scalar or 1xN vector [P.Asat]
%     .a        training amplification           [P.p]
%     .shrink   scalar multiplying the estimate  [1]     (BLMMSE)
%     .N,.K,.T  override the defaults in P
%     .sumMSE   true -> report sum-NMSE over [h_d ; b]   [false]
%
%   The DFT training design uses columns 1..N: the DC column is reserved for
%   the direct link, otherwise it is collinear with it and LS is singular.

if ~isfield(opt,'model'),  opt.model  = 'rapp';  end
if ~isfield(opt,'PT'),     opt.PT     = 1000;    end
if ~isfield(opt,'gain'),   opt.gain   = P.p;     end
if ~isfield(opt,'Asat'),   opt.Asat   = P.Asat;  end
if ~isfield(opt,'a'),      opt.a      = P.p;     end
if ~isfield(opt,'shrink'), opt.shrink = 1;       end
if ~isfield(opt,'sumMSE'), opt.sumMSE = false;   end
if isfield(opt,'N'), P.N = opt.N; end
if isfield(opt,'K'), P.K = opt.K; end
if isfield(opt,'T'), P.T = opt.T; else, P.T = P.N + 1; end

N = P.N; K = P.K; T = P.T; PT = opt.PT;
F = exp(-2j*pi*(0:T-1).' * (1:N) / T);      % T x N, unit modulus

num = 0; den = 0; nrm = 0;
for r = 1:P.nReal
    b  = sqrt(P.rho_b/2)*(randn(N,1) + 1j*randn(N,1));
    hd = sqrt(P.rho_d/2)*(randn(K,1) + 1j*randn(K,1));

    % RIS -> BS: rank-one line of sight, constant modulus per element.
    ar = exp(1j*pi*(0:K-1).' * sin(rand*pi - pi/2));
    at = exp(1j*pi*(0:N-1).' * sin(rand*pi - pi/2));
    G  = sqrt(P.rho_g) * (ar * at');                    % K x N

    v = sqrt(P.sig1_2/2)*(randn(T,N) + 1j*randn(T,N));
    z = sqrt(P.sig2_2/2)*(randn(T,K) + 1j*randn(T,K));

    u = sqrt(PT)*repmat(b.',T,1) + v;                   % element input
    switch lower(opt.model)
        case 'linear', out = opt.a * u;
        case 'rapp',   out = rapp(u, opt.a, opt.Asat, P.smooth);
        case 'saleh',  out = saleh(u, opt.a, opt.Asat);
        otherwise,     error('unknown model');
    end
    Y = (out .* F) * G.' + sqrt(PT)*repmat(hd.',T,1) + z;

    bh = zeros(N,1);
    for k = 1:K
        Phi = [ones(T,1), opt.gain * (F .* repmat(G(k,:),T,1))];
        est = (Phi \ Y(:,k)) / sqrt(PT);
        bh  = bh + est(2:end);
        if opt.sumMSE
            % per-antenna accumulation: the direct link has one coefficient
            % per antenna, so it must not be averaged across antennas first.
            num = num + abs(hd(k) - est(1))^2 + sum(abs(b - est(2:end)).^2)/K;
            den = den + abs(hd(k))^2          + sum(abs(b).^2)/K;
        end
    end
    bh = opt.shrink * bh / K;

    if ~opt.sumMSE
        num = num + sum(abs(b - bh).^2);
        den = den + sum(abs(b).^2);
    end
    nrm = nrm + norm(bh);
end
nmse      = num/den;
bhat_norm = nrm/P.nReal;
end
