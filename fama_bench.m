function P = fama_bench(which,U,tauS,cfg)
%FAMA_BENCH  Competing closed-form floors, recoded for comparison.
%
%   'edof'  Equivalent-degrees-of-freedom FAMA floor (Rayleigh, equal-power
%           emitters, perfect CSI): sum_j C(K,j)(-1)^j (1+j*tau)^{-(U-1)},
%           with K = 2*ceil(W)+1 independent modes.
%           VERIFIED to coincide with Theorem 1 at B=1, m_b -> inf, to 8
%           significant digits (see fama_selftest, test 5).
%
%   'bcm'   Block-independent Nakagami baseline in the style of the spatial
%           block-correlation literature: equal-power emitters, perfect CSI,
%           blocks treated as independent, no texture coupling.
%
%   NOTE: 'bcm' is a reconstruction of the modelling assumptions, not a
%   transcription of any particular published equation. Verify against the
%   source before submission.
if nargin < 4, cfg = fama_config(); end
if nargin < 3 || isempty(tauS), tauS = cfg.tau; end
switch lower(which)
  case 'edof'
    K = 2*ceil(cfg.W)+1;
    j = (0:K).';
    lg = gammaln(K+1)-gammaln(j+1)-gammaln(K-j+1);
    P  = sum((-1).^j .* exp(lg) .* (1+j*tauS).^(-(U-1)));
  case 'bcm'
    % Uniform severity across blocks, the assumption made when a single
    % Nakagami parameter is fitted to the whole aperture. Everything else
    % matches Corollary 2, so the curve isolates the cost of that assumption.
    mbar = mean(cfg.mb);
    [Ibar,mIb] = fama_profile(U,'hom',cfg);
    P = 1;
    for b = 1:cfg.B
        xb = mbar*Ibar(b)*tauS/(mbar*Ibar(b)*tauS + mIb(b)*cfg.Omb(b));
        P  = P*betainc(xb,mbar,mIb(b));
    end
  otherwise
    error('fama_bench:which','which must be edof or bcm');
end
end
