function [xi,se2,Sidx] = fama_lmmse(b,gp,cfg)
%FAMA_LMMSE  Captured and residual power of the block LMMSE estimator.
%   Uses the PHYSICAL correlation of eq. (5) restricted to block b, and a
%   uniformly spaced probing subset S_b. Returns the two second-order
%   quantities the DP-PSINR metric consumes, eq. (11).
if nargin<3, cfg = fama_config(); end
persistent J
if isempty(J), J = fama_aperture(cfg); end
s0 = sum(cfg.Lb(1:b-1)); idx = s0+(1:cfg.Lb(b));
Rb = cfg.Omb(b)*real(J(idx,idx));
nS = max(2,round(cfg.Lb(b)/3));
Sidx = round(linspace(1,cfg.Lb(b),nS));
Rss = Rb(Sidx,Sidx); Rcs = Rb(:,Sidx);
M   = gp*Rss + eye(nS);                      % normalized by sigma_n^2
Cap = gp*(Rcs/M)*Rcs.';                      % captured covariance
xi  = trace(Cap)/cfg.Lb(b);
se2 = trace(Rb)/cfg.Lb(b) - xi;
end
