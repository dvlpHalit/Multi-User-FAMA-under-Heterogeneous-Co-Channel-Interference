function fig5_diversity()
%FIG5  Floor diversity: slope -> d_inf = sum_b m_b, independent of L_b.
%      Uses Corollary 2, which carries no alternating sum and is therefore
%      the only expression evaluable in double precision as tau_S -> 0.
cfg = fama_config(); ts = logspace(-1,-6,26);
h = figure; hold on; lg = {};
for U = [3 5 8]
    P = arrayfun(@(t) fama_floor('coarse',U,'hom',t,cfg), ts);
    loglog(ts,P,'-','LineWidth',1.3);
    p = polyfit(log10(ts),log10(P),1);
    lg{end+1} = sprintf('$U=%d$, slope $=%.2f$',U,p(1));
end
loglog(ts,1e-2*(ts/ts(1)).^cfg.dinf,'k--','LineWidth',1.0);
lg{end+1} = sprintf('$\\tau_{\\mathrm{S}}^{d_{\\infty}}$, $d_{\\infty}=%.1f$',cfg.dinf);
set(gca,'XScale','log','YScale','log');
fama_style(h,'$\tau_{\mathrm{S}}$','$P_{\infty}$',lg,'southeast');
fama_print(h,'fig5_diversity');
end
