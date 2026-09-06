function cfg = fama_config()
%FAMA_CONFIG  Single source of truth for every simulation parameter.
%   Deliberately distinct from the thesis configuration so that no figure in
%   the paper duplicates a thesis figure.
%
%   thesis : N=64, W=2.5, kappaA=5, B=4, L_b=16, m_b=2.5 (uniform), tau=3 dB
%   paper  : N=100, W=3.0, kappaA=8, B=4, L_b=25, m_b heterogeneous, tau=0 dB

cfg.N       = 100;                 % candidate ports
cfg.W       = 3.0;                 % aperture, wavelengths
cfg.kappaA  = 8.0;                 % von Mises concentration
cfg.dx      = cfg.W/(cfg.N-1);     % port spacing, wavelengths
cfg.mu      = 0.0;                 % mean AoA (irrelevant: Proposition 1)

cfg.B       = 4;                   % blocks  (set from aperture r_eff, see fama_aperture)
cfg.Lb      = [25 25 25 25];       % ports per block
cfg.mb      = [3.0 2.4 1.8 1.4];   % per-block severity -> d_inf = 8.6
cfg.Omb     = [1.0 1.0 1.0 1.0];   % per-block mean power

cfg.mI      = 2.5;                 % emitter severity
cfg.SIRdB   = 10;                  % per-emitter SIR, dB
cfg.taudB   = 0;                   % outage threshold, dB
cfg.eta     = 3.5;                 % path-loss exponent (heterogeneous profile)

cfg.Uvec    = 2:8;                 % emitter counts swept
cfg.nGL      = 64;                 % generalized Gauss-Laguerre nodes
cfg.nPanel   = 46;                 % panels in the adaptive texture grid
cfg.nNode    = 20;                 % Gauss-Legendre nodes per panel
cfg.seed     = 20260901;

cfg.dinf    = sum(cfg.mb);
cfg.tau     = 10^(cfg.taudB/10);
end
