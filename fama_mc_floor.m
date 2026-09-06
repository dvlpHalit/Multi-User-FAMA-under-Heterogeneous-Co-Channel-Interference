function [P,ci] = fama_mc_floor(model,U,mode,tauS,ntr,cfg)
%FAMA_MC_FLOOR  Monte Carlo outage floor.
%   model 'block' : the model the analysis assumes. Block-common desired
%                   texture, block-common interference, ports conditionally
%                   independent. Validates Theorem 1.
%   model 'full'  : the physical field. Ports correlated through the von
%                   Mises matrix, emitter textures APERTURE-common. This is
%                   P_inf^fc, the ground truth of Section V.
if nargin < 6, cfg = fama_config(); end
if nargin < 5 || isempty(ntr),  ntr  = 2e5; end
if nargin < 4 || isempty(tauS), tauS = cfg.tau; end
rand('twister',cfg.seed); randn('state',cfg.seed);

[Ibar,mIb,~,Omv] = fama_profile(U,mode,cfg);
B = cfg.B; Lb = cfg.Lb; mb = cfg.mb; Omb = cfg.Omb;
chunk = 2e4; done = 0; nout = 0;

if strcmpi(model,'full')
    [J,~,~] = fama_aperture(cfg);
    idx = cell(1,B); s = 0;
    for b = 1:B, idx{b} = s+(1:Lb(b)); s = s+Lb(b); end
    Lc  = chol(J + 1e-10*eye(cfg.N));           % upper, J = Lc'*Lc
end

while done < ntr
    n = min(chunk,ntr-done);
    if strcmpi(model,'block')
        fail = true(n,1);
        for b = 1:B
            t  = fama_gamrnd(mb(b),1/mb(b),[n 1]);
            I  = fama_gamrnd(mIb(b),Ibar(b)/mIb(b),[n 1]);
            e  = -log(rand(n,Lb(b)));                       % Exp(1)
            S  = Omb(b)*(t*ones(1,Lb(b))).*e;
            fail = fail & all(S < tauS*(I*ones(1,Lb(b))),2);
        end
    else
        G = (randn(n,cfg.N)+1i*randn(n,cfg.N))/sqrt(2);
        Sd = zeros(n,cfg.N);
        for b = 1:B
            Jb = J(idx{b},idx{b});
            Cb = chol(Jb + 1e-10*eye(Lb(b)));
            gb = G(:,idx{b})*Cb;
            t  = fama_gamrnd(mb(b),1/mb(b),[n 1]);
            Sd(:,idx{b}) = Omb(b)*(t*ones(1,Lb(b))).*abs(gb).^2;
        end
        Itot = zeros(n,cfg.N);
        for v = 1:numel(Omv)
            Gv = ((randn(n,cfg.N)+1i*randn(n,cfg.N))/sqrt(2))*Lc;
            tv = fama_gamrnd(cfg.mI,1/cfg.mI,[n 1]);        % APERTURE-common
            Itot = Itot + Omv(v)*(tv*ones(1,cfg.N)).*abs(Gv).^2;
        end
        fail = all(Sd < tauS*Itot,2);
    end
    nout = nout + sum(fail); done = done + n;
end
P  = nout/ntr;
ci = 1.96*sqrt(max(nout,1))/ntr;
end
