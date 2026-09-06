function [Ibar,mIb,Ueff,Omv] = fama_profile(U,mode,cfg)
%FAMA_PROFILE  Emitter power profile -> (Ibar_b, m_I,b, U_eff,b) via Lemma 1.
%   mode 'hom' : equal-power emitters                    -> U_eff = U-1
%   mode 'het' : path-loss spread over d in [1,4]        -> U_eff < U-1
%   Both modes are normalized to the SAME total mean interference Ibar, so
%   Fig. 6 isolates the effect of the profile shape at fixed Ibar.
if nargin < 3, cfg = fama_config(); end
V   = U-1;
Om0 = 10^(-cfg.SIRdB/10);                 % per-emitter power in the 'hom' case
if V < 1, Ibar = zeros(1,cfg.B); mIb = ones(1,cfg.B); Ueff = 0; Omv = []; return; end
switch lower(mode)
  case 'hom'
    Omv = Om0*ones(1,V);
  case 'het'
    d   = linspace(1,4,V);
    Omv = d.^(-cfg.eta);
    Omv = Omv*(V*Om0)/sum(Omv);           % same total power as 'hom'
  otherwise
    error('fama_profile:mode','mode must be hom or het');
end
Ueff = sum(Omv)^2/sum(Omv.^2);            % Definition 1
Ibar = sum(Omv)*ones(1,cfg.B);            % flat in-aperture profile
mIb  = (sum(Omv)^2/sum(Omv.^2/cfg.mI))*ones(1,cfg.B);   % Lemma 1
end
