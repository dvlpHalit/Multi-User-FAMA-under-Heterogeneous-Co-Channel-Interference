function fama_selftest()
%FAMA_SELFTEST  Verifies every numerical routine before any figure is made.
cfg = fama_config(); np = 0; nf = 0;
fprintf('\n=== FAMA numerical self-test ===\n');

  function chk(name,got,want,tol)
    e = abs(got/want-1);
    if e <= tol, fprintf('  PASS  %-46s err=%.1e\n',name,e); np=np+1;
    else,        fprintf('  FAIL  %-46s err=%.1e\n',name,e); nf=nf+1; end
  end

% 1 quadrature rules integrate their own moments
[x,w] = fama_gauss('legendre01',40);
chk('Gauss-Legendre(0,1): int x^3 = 1/4', sum(w.*x.^3), 0.25, 1e-12);
[x,w] = fama_gauss('laguerre',40,1.5);
chk('gen. Gauss-Laguerre: E[X] for Gamma(2.5,1)', sum(w.*x), 2.5, 1e-12);
chk('gen. Gauss-Laguerre: E[X^2]', sum(w.*x.^2), 2.5*3.5, 1e-11);

% 2 aperture degrees of freedom set the block count
[~,reff,B] = fama_aperture(cfg);
chk('aperture r_eff gives the configured B', B, cfg.B, 1e-12);
fprintf('        r_eff = %.3f  ->  B = %d,  L_b = %d\n',reff,B,cfg.Lb(1));

% 3 block floor against the CLOSED FORM of Theorem 1 (Tricomi), small L so
%   that the alternating sum still has usable precision
L = 6; m = 2.0; mI = 3.0; Ib = 0.4; Om = 1.0; tauS = 0.5;
ab = tauS*Ib/(mI*Om);  ser = 1;
for l = 1:L
    lg = gammaln(L+1)-gammaln(l+1)-gammaln(L-l+1);
    T  = exp(gammaln(m+mI)-gammaln(m)+m*log(m)+m*log(l*ab)) ...
         * fama_tricomi(m+mI,m+1,m*l*ab);
    ser = ser + (-1)^l*exp(lg)*T;
end
[u2,w2] = fama_gauss('laguerre',200,mI-1);
qd = sum(w2.*fama_blockfloor(m,L,(tauS*(Ib/mI)*u2/Om).',cfg).');
chk('Theorem 1 closed form vs quadrature (L=6)', qd, ser, 1e-6);

% 4 consistency reductions of Section IV-D
chk('(S1) floor vanishes as Ibar -> 0', ...
    fama_floor('coarse',1,'hom',cfg.tau,cfg)+1, 1, 1e-12);
Pc = fama_floor('coarse',4,'hom',cfg.tau,cfg);
Pt = fama_floor('texture',4,'hom',cfg.tau,cfg);
if Pt >= Pc, fprintf('  PASS  %-46s\n','texture floor >= coarse floor'); np=np+1;
else,        fprintf('  FAIL  %-46s\n','texture floor >= coarse floor'); nf=nf+1; end

% 5 the EDoF benchmark is Theorem 1 at B=1, m_b -> inf
c1 = cfg; c1.B = 1; c1.Lb = 2*ceil(cfg.W)+1; c1.mb = 1e7; c1.Omb = 1;
for U = [3 5]
    K = 2*ceil(cfg.W)+1; j = (0:K).';
    lg = gammaln(K+1)-gammaln(j+1)-gammaln(K-j+1);
    wu = sum((-1).^j.*exp(lg).*(1+j*cfg.tau).^(-(U-1)));
    c1.mI = 1;                         % mIb = (U-1)*mI = U-1  (Rayleigh emitters)
    c1.SIRdB = 0;                      % per-emitter power 1 -> Ibar = U-1
    ours = fama_floor('fine',U,'hom',cfg.tau,c1);
    chk(sprintf('EDoF floor == Thm 1 (B=1, m->inf), U=%d',U), ours, wu, 5e-4);
end

% 6 floor diversity order d_inf = sum_b m_b (cancellation-free branch)
ts = logspace(-4,-6,5); v = arrayfun(@(t) fama_floor('coarse',4,'hom',t,cfg), ts);
pf = polyfit(log10(ts),log10(v),1);
chk('floor diversity slope = sum_b m_b', pf(1), cfg.dinf, 2e-3);

% 7 Tricomi cross-check
chk('Tricomi U(5,3.5,2) integral identity', ...
    fama_tricomi(5,3.5,2), 2.08596257e-03, 1e-6);

% 8 gamma sampler moments
g = fama_gamrnd(2.5,1/2.5,[4e5 1]);
chk('Marsaglia-Tsang mean',     mean(g), 1.0,     2e-2);
chk('Marsaglia-Tsang variance', var(g),  1/2.5,   5e-2);

fprintf('=== %d passed, %d failed ===\n\n',np,nf);
if nf>0, error('fama_selftest:failed','%d test(s) failed',nf); end
end
