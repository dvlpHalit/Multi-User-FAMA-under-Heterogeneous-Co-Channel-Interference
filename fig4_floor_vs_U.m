function fig4_floor_vs_U()
%FIG4  Floor versus emitter count: the three expressions, the correlated
%      field, and the two external baselines.  Panel (b): amplification
%      budget against the FGM structural ceiling 2^{B-1}.
cfg = fama_config(); U = cfg.Uvec; nU = numel(U);
Pc=zeros(1,nU); Pf=Pc; Pt=Pc; Pm=Pc; ci=Pc; Pe=Pc; Pb=Pc;
for k = 1:nU
    Pc(k) = fama_floor('coarse' ,U(k),'hom',cfg.tau,cfg);
    Pf(k) = fama_floor('fine'   ,U(k),'hom',cfg.tau,cfg);
    Pt(k) = fama_floor('texture',U(k),'hom',cfg.tau,cfg);
    nt = 2e5; if U(k)<=3, nt = 1.2e6; end   % deep floors need more events
    [Pm(k),ci(k)] = fama_mc_floor('full',U(k),'hom',cfg.tau,nt,cfg);
    Pe(k) = fama_bench('edof',U(k),cfg.tau,cfg);
    Pb(k) = fama_bench('bcm' ,U(k),cfg.tau,cfg);
    fprintf('  U=%d done\n',U(k));
end
h = figure;
subplot(1,2,1); set(gca,'YScale','log'); hold on;
plot(U,Pf,'v-','LineWidth',1.1); plot(U,Pc,'s-','LineWidth',1.1);
plot(U,Pt,'d-','LineWidth',1.4);
% clamp the lower whisker: on a log axis a negative bound is silently dropped
errorbar(U,Pm,min(ci,0.9*Pm),ci,'ko','MarkerSize',4);
plot(U,Pe,'--','LineWidth',1.0); plot(U,Pb,':','LineWidth',1.0);
fama_style(h,'number of co-channel links $U$','$P_{\infty}$', ...
  {'Theorem 1 (fine)','Corollary 2 (compact)','Theorem 4 (texture)', ...
   'correlated field (MC)','EDoF-FAMA','BCM baseline'},'southeast');
subplot(1,2,2); set(gca,'YScale','log'); hold on;
plot(U,Pm./Pc,'ko-','LineWidth',1.2);
plot(U,Pt./Pc,'d-','LineWidth',1.4);
plot(U,2^(cfg.B-1)*ones(1,nU),'r--','LineWidth',1.2);
fama_style(h,'number of co-channel links $U$','amplification over $P_{\infty}$', ...
  {'required, $\mathcal{A}$','supplied by Theorem 4','FGM ceiling $2^{B-1}$'},'northeast');
set(h,'Units','inches','Position',[1 1 7.16 2.75]);   % double column
fama_print(h,'fig4_floor_vs_U');
save('fig4_data.mat','U','Pc','Pf','Pt','Pm','ci','Pe','Pb');
fprintf('\n  Fig.4 data (paste this back):\n');
fprintf('   U  coarse      fine        texture     MCfull      +/-\n');
for k=1:nU, fprintf('   %d  %.4e  %.4e  %.4e  %.4e  %.1e\n',U(k),Pc(k),Pf(k),Pt(k),Pm(k),ci(k)); end
end
