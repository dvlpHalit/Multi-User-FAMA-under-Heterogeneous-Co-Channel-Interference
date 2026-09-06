function [J,reff,B] = fama_aperture(cfg)
%FAMA_APERTURE  Physical von Mises / Abdi-Kaveh port correlation, eq. (5),
%   its effective rank, eq. (8), and the implied block count.
%   Simulations use this PHYSICAL matrix rather than the banded surrogate, so
%   the analysis is validated against the model it approximates.
if nargin < 1, cfg = fama_config(); end
i  = (0:cfg.N-1).';
D  = 2*pi*cfg.dx*abs(i - i.');
arg = sqrt(cfg.kappaA^2 - D.^2 + 1i*2*cfg.kappaA*D*sin(cfg.mu));
J   = besseli(0,arg)/besseli(0,cfg.kappaA);
J   = (J + J')/2;                                  % enforce Hermitian
reff = real(trace(J))^2/real(trace(J*J'));
B    = round(reff);
end
