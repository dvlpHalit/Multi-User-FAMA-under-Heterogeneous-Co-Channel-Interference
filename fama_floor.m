function P = fama_floor(which,U,mode,tauS,cfg)
%FAMA_FLOOR  Interference-limited outage floor, three expressions.
%   which = 'coarse'   Corollary 2, eq. (20)     prod_b I_{x_b}(m_b,m_I,b)
%         = 'fine'     Theorem 1,   eq. (18)     cancellation-free evaluation
%         = 'texture'  Theorem 4,   eq. (31)     aperture-common texture
if nargin < 5, cfg = fama_config(); end
if nargin < 4 || isempty(tauS), tauS = cfg.tau; end
[Ibar,mIb] = fama_profile(U,mode,cfg);
mb = cfg.mb; Omb = cfg.Omb; Lb = cfg.Lb; B = cfg.B;

switch lower(which)
  case 'coarse'
    P = 1;
    for b = 1:B
        xb = mb(b)*Ibar(b)*tauS/(mb(b)*Ibar(b)*tauS + mIb(b)*Omb(b));
        P  = P*betainc(xb,mb(b),mIb(b));
    end

  case 'fine'
    P = 1;
    for b = 1:B
        [u,w] = fama_gauss('laguerre',cfg.nGL,mIb(b)-1);   % u ~ Gamma(mIb,1)
        Ik    = (Ibar(b)/mIb(b))*u;                        % emitter power draws
        ck    = tauS*Ik/Omb(b);
        Fk    = fama_blockfloor(mb(b),Lb(b),ck,cfg);
        P     = P*sum(w(:).*Fk(:));
    end

  case 'texture'
    mUps = (U-1)*cfg.mI;                                   % eq. (30)
    cb   = (mUps+1)/(cfg.mI+1);
    [u,w] = fama_gauss('laguerre',cfg.nGL,mUps-1);
    ups   = u/mUps;                                        % ups ~ Gamma(mUps,1/mUps)
    inner = ones(size(ups));
    for b = 1:B
        xb    = mb(b)*ups*Ibar(b)*tauS ./ (mb(b)*ups*Ibar(b)*tauS + cb*Omb(b));
        inner = inner.*betainc(xb,mb(b),cb);
    end
    P = sum(w(:).*inner(:));

  otherwise
    error('fama_floor:which','which must be coarse, fine or texture');
end
end
