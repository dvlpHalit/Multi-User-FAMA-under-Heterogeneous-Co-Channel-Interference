function fig3_outage_snr()
%FIG3  Outage versus pilot SNR: emergence of the interference-limited floor.
%      Selection uses the deterministic-equivalent surrogate (16); scoring
%      uses the actual DP-PSINR (15), as required by Remark 1.
cfg = fama_config(); rand('twister',cfg.seed); randn('state',cfg.seed);
gpdB = 0:5:45; gp = 10.^(gpdB/10); Uset = [1 2 8]; ntr = 4e4;
J = fama_aperture(cfg); B = cfg.B; Lb = cfg.Lb;
idx = cell(1,B); s0 = 0;
for b=1:B, idx{b} = s0+(1:Lb(b)); s0 = s0+Lb(b); end
Cb = cell(1,B);
for b=1:B, Cb{b} = chol(real(J(idx{b},idx{b})) + 1e-10*eye(Lb(b))); end
Lc = chol(real(J) + 1e-10*eye(cfg.N));
P = zeros(numel(Uset),numel(gp));
for iu = 1:numel(Uset)
    U = Uset(iu); [~,~,~,Omv] = fama_profile(U,'hom',cfg);
    for ig = 1:numel(gp)
        g = gp(ig); nout = 0;
        xi = zeros(1,B); se2 = zeros(1,B);
        for b=1:B, [xi(b),se2(b)] = fama_lmmse(b,g,cfg); end
        Psi = zeros(ntr,cfg.N); Psih = zeros(ntr,cfg.N);
        Itot = zeros(ntr,cfg.N);
        for v = 1:numel(Omv)
            Gv = ((randn(ntr,cfg.N)+1i*randn(ntr,cfg.N))/sqrt(2))*Lc;
            tv = fama_gamrnd(cfg.mI,1/cfg.mI,[ntr 1]);
            Itot = Itot + Omv(v)*(tv*ones(1,cfg.N)).*abs(Gv).^2;
        end
        Ibar = sum(Omv);
        for b = 1:B
            t  = fama_gamrnd(cfg.mb(b),1/cfg.mb(b),[ntr 1]);
            gb = ((randn(ntr,Lb(b))+1i*randn(ntr,Lb(b)))/sqrt(2))*Cb{b};
            ph = xi(b)*(t*ones(1,Lb(b))).*abs(gb).^2;      % |ghat|^2, mean xi_b
            den  = g*se2(b) + g*Itot(:,idx{b}) + 1;
            denh = g*se2(b) + g*Ibar + 1;
            Psi (:,idx{b}) = g*ph./den;
            Psih(:,idx{b}) = g*ph./denh;
        end
        [~,sel] = max(Psih,[],2);
        lin = sub2ind(size(Psi),(1:ntr).',sel);
        P(iu,ig) = mean(Psi(lin) < cfg.tau);
    end
    fprintf('  U=%d done\n',U);
end
h = figure; hold on; lg = {}; C = lines(numel(Uset));
for iu = 1:numel(Uset)
    semilogy(gpdB,max(P(iu,:),1e-6),'o-','Color',C(iu,:),'LineWidth',1.2);
    lg{end+1} = sprintf('$U=%d$',Uset(iu));
    if Uset(iu) > 1
        Pf = fama_floor('texture',Uset(iu),'hom',cfg.tau,cfg);
        semilogy(gpdB,Pf*ones(size(gpdB)),'--','Color',C(iu,:),'LineWidth',0.9);
        lg{end+1} = sprintf('floor, Thm 4, $U=%d$',Uset(iu));
    end
end
set(gca,'YScale','log');
fama_style(h,'pilot SNR $\gamma_p$ (dB)','$P_{\mathrm{out}}(\tau)$',lg,'northeast');
fama_print(h,'fig3_outage_snr');
save('fig3_data.mat','gpdB','P','Uset');
fprintf('\n  Fig.3 sweep (paste this back):\n   gp(dB) ');
for iu=1:numel(Uset), fprintf('     U=%d',Uset(iu)); end
fprintf('\n');
for ig=1:numel(gpdB)
    fprintf('   %5d ',gpdB(ig)); fprintf('%9.2e',P(:,ig)); fprintf('\n');
end
fprintf('   floor(Thm4) ');
for iu=1:numel(Uset)
    if Uset(iu)>1, fprintf('%9.2e',fama_floor('texture',Uset(iu),'hom',cfg.tau,cfg));
    else, fprintf('%9s','-'); end
end
fprintf('\n');
end
