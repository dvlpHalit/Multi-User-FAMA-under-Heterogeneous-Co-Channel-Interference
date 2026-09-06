function U = fama_tricomi(a,c,z,n)
%FAMA_TRICOMI  Tricomi confluent hypergeometric function of the second kind.
%   Base MATLAB has no hyperu (Symbolic Toolbox provides kummerU). Two
%   branches of the integral representation
%       U(a,c,z) = 1/Gamma(a) int_0^inf e^{-zt} t^{a-1} (1+t)^{c-a-1} dt
%   are used: substitution t=u/z with generalized Gauss-Laguerre for z>=1,
%   and a finite-interval map with the endpoint singularity removed for z<1.
%
%   Used only to cross-check the closed form of Theorem 1 in fama_selftest;
%   the figures use fama_blockfloor, which is cancellation-free.
%   Accuracy is ~1e-9 for z>=0.02 and degrades for very small z with a<3.
if nargin < 4, n = 200; end
if z >= 1
    [x,w] = fama_gauss('laguerre',n,a-1);
    U = z^(-a)*sum(w.*(1+x/z).^(c-a-1));
else
    [v,w] = fama_gauss('legendre01',n);
    q = 1 - v.^(1/a);
    U = sum(w.*exp(-z*(1-q)./q - c*log(q)))/(a*gamma(a));
end
end
