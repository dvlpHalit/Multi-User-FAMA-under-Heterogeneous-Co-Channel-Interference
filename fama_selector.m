function out = fama_selector(U,cbar,gpdB,ntr,cfg)
%FAMA_SELECTOR  Cost of ignoring the interference profile in port selection.
%   Compares three selectors, all SCORED on the actual DP-PSINR of eq. (15):
%     aware    : eq. (16), uses the per-block mean interference Ibar_b
%     agnostic : maximizes the estimated signal term alone
%     oracle   : instantaneous best port (upper bound, needs full CSI)
%   cbar is the normalized in-aperture interference profile of eq. (10);
%   a flat cbar makes the two selectors identical, which is Proposition 3
%   with nothing left to exploit.
if nargin<5, cfg=fama_config(); end
rand('twister',cfg.seed); randn('state',cfg.seed);
gp = 10^(gpdB/10); B=cfg.B; Lb=cfg.Lb;
J = fama_aperture(cfg); s0=0; idx=cell(1,B); Cb=cell(1,B);
for b=1:B, idx{b}=s0+(1:Lb(b)); s0=s0+Lb(b); Cb{b}=chol(real(J(idx{b},idx{b}))+1e-10*eye(Lb(b))); end
Lc = chol(real(J)+1e-10*eye(cfg.N));
[~,~,~,Omv] = fama_profile(U,'hom',cfg);
xi=zeros(1,B); se2=zeros(1,B);
for b=1:B, [xi(b),se2(b)]=fama_lmmse(b,gp,cfg); end
cvec = zeros(1,cfg.N);
for b=1:B, cvec(idx{b}) = cbar(b); end
Itot = zeros(ntr,cfg.N);
for v=1:numel(Omv)
    Gv = ((randn(ntr,cfg.N)+1i*randn(ntr,cfg.N))/sqrt(2))*Lc;
    tv = fama_gamrnd(cfg.mI,1/cfg.mI,[ntr 1]);
    Itot = Itot + Omv(v)*(tv*ones(1,cfg.N)).*abs(Gv).^2.*(ones(ntr,1)*cvec);
end
Ibar_b = sum(Omv)*cvec;                     % deterministic equivalent, eq. (10)
Psi=zeros(ntr,cfg.N); Sig=zeros(ntr,cfg.N); Pag=zeros(ntr,cfg.N);
for b=1:B
    t  = fama_gamrnd(cfg.mb(b),1/cfg.mb(b),[ntr 1]);
    gb = ((randn(ntr,Lb(b))+1i*randn(ntr,Lb(b)))/sqrt(2))*Cb{b};
    ph = xi(b)*(t*ones(1,Lb(b))).*abs(gb).^2;
    Sig(:,idx{b}) = ph;
    Psi(:,idx{b}) = gp*ph./(gp*se2(b)+gp*Itot(:,idx{b})+1);
    Pag(:,idx{b}) = gp*ph./(gp*se2(b)+1);                 % ignores interference
end
Paw = gp*Sig./(gp*(ones(ntr,1)*(se2(ceil((1:cfg.N)/Lb(1))))) ...
      + gp*(ones(ntr,1)*Ibar_b) + 1);                     % eq. (16)
pick = @(M) sub2ind(size(Psi),(1:ntr).',max_idx(M));
out.aware    = mean(Psi(pick(Paw)) < cfg.tau);
out.agnostic = mean(Psi(pick(Pag)) < cfg.tau);
out.oracle   = mean(max(Psi,[],2)  < cfg.tau);
end
function k = max_idx(M), [~,k] = max(M,[],2); end
