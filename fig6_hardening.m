function fig6_hardening()
%FIG6  Interference hardening: homogeneous vs heterogeneous profiles at the
%      SAME mean interference power, and the two-directional threshold effect.
cfg = fama_config(); U = cfg.Uvec;
Ph=zeros(size(U)); Pt=Ph; Ue=Ph;
for k=1:numel(U)
    Ph(k) = fama_floor('coarse',U(k),'hom',cfg.tau,cfg);
    Pt(k) = fama_floor('coarse',U(k),'het',cfg.tau,cfg);
    [~,~,Ue(k)] = fama_profile(U(k),'het',cfg);
end
h = figure;
subplot(1,2,1); hold on;
semilogy(U,Ph,'s-','LineWidth',1.3); semilogy(U,Pt,'d-','LineWidth',1.3);
set(gca,'YScale','log');
fama_style(h,'number of co-channel links $U$','$P_{\infty}$', ...
  {'homogeneous, $U_{\mathrm{eff}}=U-1$','heterogeneous, same $\bar{I}_b$'},'southeast');
subplot(1,2,2); hold on;
plot(U,U-1,'s-','LineWidth',1.3); plot(U,Ue,'d-','LineWidth',1.3);
fama_style(h,'number of co-channel links $U$','$U_{\mathrm{eff},b}$', ...
  {'homogeneous','heterogeneous'},'northwest');
set(h,'Units','inches','Position',[1 1 7.16 2.75]);
fama_print(h,'fig6_hardening');
end
