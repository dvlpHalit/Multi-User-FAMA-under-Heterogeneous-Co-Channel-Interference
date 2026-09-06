function F = fama_blockfloor(m,L,c,cfg)
%FAMA_BLOCKFLOOR  E_t[(1-exp(-c/t))^L] for t ~ Gamma(m,1/m).
%   Per-block interference-limited outage factor of Theorem 1, evaluated
%   WITHOUT the alternating binomial sum. That sum loses up to 20 significant
%   digits at tau_S = 1e-5 and is unusable in double precision below
%   tau_S ~ 1e-3; this routine is cancellation-free.
%
%   Method: substitute t = F_t^{-1}(p) and integrate over p in (0,1) on a
%   panel grid anchored at the half-power texture t*, where the block outage
%   factor equals 1/2, so resolution follows the floor depth.
%
%   Validated to ~1e-15 relative error for m in [1.2,3.0], c in [1e-6,1e-1]
%   (fama_selftest, test 3). For m > 1e4 the texture is deterministic and the
%   closed limit (1-exp(-c))^L is returned instead, which also avoids the
%   very slow gammaincinv evaluations at large shape.
if nargin < 4, cfg = fama_config(); end
c  = c(:).';
F  = zeros(size(c));
if m > 1e4                                  % t -> 1 in probability
    F = exp(L*log(-expm1(-c)));
    F(c<=0) = 0;
    return
end
[x0,w0] = fama_gauss('legendre01',cfg.nNode);
x0 = x0(:).'; w0 = w0(:).';                 % row vectors for the outer product
for q = 1:numel(c)
    cq = c(q);
    if cq <= 0, continue; end
    tstar = cq/(-log1p(-2^(-1/L)));         % (1-exp(-cq/t*))^L = 1/2
    pstar = gammainc(m*tstar,m);
    plo   = max(min(pstar*1e-14,1e-2),1e-300);
    edges = [0, logspace(log10(plo),0,cfg.nPanel)];
    lo = edges(1:end-1).';  hi = edges(2:end).';
    X  = reshape((lo*ones(1,cfg.nNode) + (hi-lo)*x0).',[],1);
    Wt = reshape(((hi-lo)*w0).',[],1);
    t  = gammaincinv(X,m)/m;
    g  = exp(L*log(-expm1(-cq./t)));
    g(~isfinite(g)) = 0;
    F(q) = sum(Wt.*g);
end
end
