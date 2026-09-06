function fig2_sircdf()
%FIG2  Validation of the Welch-Satterthwaite reduction (Lemma 1 -> Lemma 2).
%      Block SIR CDF, analytic F law vs Monte Carlo, homogeneous and
%      heterogeneous emitter profiles at the SAME mean interference power.
cfg = fama_config(); rand('twister',cfg.seed); randn('state',cfg.seed);
b = 1; U = 6; n = 4e5;
z = logspace(-1.5,1.5,60).';
h = figure; hold on; C = lines(4); k = 0; lg = {};
for mode = {'hom','het'}
    k = k+1;
    [Ibar,mIb,Ueff,Omv] = fama_profile(U,mode{1},cfg);
    xz  = cfg.mb(b)*Ibar(b)*z ./ (cfg.mb(b)*Ibar(b)*z + mIb(b)*cfg.Omb(b));
    Fan = betainc(xz,cfg.mb(b),mIb(b));                        % eq. (13)
    S = fama_gamrnd(cfg.mb(b),cfg.Omb(b)/cfg.mb(b),[n 1]);
    I = zeros(n,1);
    for v = 1:numel(Omv)
        I = I + fama_gamrnd(cfg.mI,Omv(v)/cfg.mI,[n 1]);       % exact sum
    end
    Zs  = sort(S./I);
    Fmc = interp1(Zs,(1:n).'/n,z,'linear','extrap');
    semilogx(z,Fan,'-','Color',C(k,:),'LineWidth',1.3);
    semilogx(z(1:4:end),Fmc(1:4:end),'o','Color',C(k,:),'MarkerSize',4);
    lg{end+1} = sprintf('%s, $U_{\\mathrm{eff}}=%.2f$ (analysis)', ...
                        upper(mode{1}),Ueff);
    lg{end+1} = sprintf('%s (simulation)',upper(mode{1}));
end
set(gca,'XScale','log');
fama_style(h,'block SIR $z$','$F_{Z_b}(z)$',lg,'northwest');
fama_print(h,'fig2_sircdf');
end
