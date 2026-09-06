function g = fama_gamrnd(shape,scale,sz)
%FAMA_GAMRND  Gamma samples by Marsaglia-Tsang. No Statistics Toolbox.
%   Handles shape < 1 by the boost g = g_{a+1} * U^{1/a}.
if numel(sz)==1, sz = [sz 1]; end
n = prod(sz);
if shape < 1
    g = fama_gamrnd(shape+1,1,[n 1]) .* rand(n,1).^(1/shape);
else
    d = shape - 1/3;  cc = 1/sqrt(9*d);
    g = zeros(n,1);   k = 0;
    while k < n
        m  = n-k;
        x  = randn(m,1);
        v  = (1+cc*x).^3;
        u  = rand(m,1);
        pos = v > 0;                       % log(v) is complex otherwise
        ok  = false(m,1);
        vp  = v(pos); xp = x(pos); up = u(pos);
        ok(pos) = log(up) < 0.5*xp.^2 + d - d*vp + d*log(vp);
        nk = sum(ok);
        if nk>0, g(k+1:k+nk) = d*v(ok); k = k+nk; end
    end
end
g = reshape(g*scale,sz);
end
