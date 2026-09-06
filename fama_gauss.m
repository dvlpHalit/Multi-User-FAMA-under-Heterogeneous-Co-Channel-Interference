function [x,w] = fama_gauss(kind,n,alpha)
%FAMA_GAUSS  Gauss rules by Golub-Welsch. No toolbox required.
%   [x,w] = fama_gauss('legendre01',n)        nodes/weights on (0,1), sum(w)=1
%   [x,w] = fama_gauss('laguerre',n,alpha)    weight x^alpha e^{-x} on (0,inf),
%                                             weights NORMALIZED to sum 1, so
%                                             sum(w.*f(x)) = E[f(X)], X~Gamma(alpha+1,1)
switch lower(kind)
  case 'legendre01'
    k  = (1:n-1)';
    b  = k./sqrt(4*k.^2-1);
    T  = diag(b,1) + diag(b,-1);
    [V,D] = eig(T);
    [x,ix] = sort(diag(D));
    w  = (V(1,ix).^2).';                 % sum = 1 on (-1,1) with weight 1/2
    x  = 0.5*(x+1);                      % map to (0,1); weights unchanged
  case 'laguerre'
    k  = (0:n-1)';
    d  = 2*k + alpha + 1;
    e  = sqrt((1:n-1)'.*((1:n-1)'+alpha));
    T  = diag(d) + diag(e,1) + diag(e,-1);
    [V,D] = eig(T);
    [x,ix] = sort(diag(D));
    w  = (V(1,ix).^2).';                 % already sums to 1
  otherwise
    error('fama_gauss:kind','unknown rule "%s"',kind);
end
w = w(:); x = x(:);
w = w/sum(w);                            % guard against round-off drift
end
