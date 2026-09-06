function collect_results_core(quick)
%COLLECT_RESULTS_CORE  Prints every number quoted in Section VII and VI.
%   quick = true  -> few Monte Carlo trials, for a smoke test only.
if nargin<1, quick=false; end
cfg = fama_config();
[~,reff,B] = fama_aperture(cfg);

fprintf('\n--- Table II: setup ---\n');
fprintf('N=%d  W=%.1f lambda  kappa_A=%.1f  dx/lambda=%.4f\n',cfg.N,cfg.W,cfg.kappaA,cfg.dx);
fprintf('r_eff = %.3f  ->  B = %d,  L_b = %d,  d_inf = %.2f\n',reff,B,cfg.Lb(1),cfg.dinf);
fprintf('m_b = [%s]   m_I = %.1f   SIR = %d dB   tau = %d dB\n', ...
        num2str(cfg.mb),cfg.mI,cfg.SIRdB,cfg.taudB);
fprintf('\nLMMSE split (Sec. VII-C), |S_b| = %d of %d:\n',max(2,round(cfg.Lb(1)/3)),cfg.Lb(1));
for g=[10 20 30]
    fprintf('  gamma_p=%2d dB: ',g);
    for b=1:cfg.B
        [xi,se2]=fama_lmmse(b,10^(g/10),cfg);
        fprintf('b%d xi=%.4f se2=%.2e  ',b,xi,se2);
    end
    fprintf('\n');
end

fprintf('\n--- Fig. 4(a) and Table III: floors vs U ---\n');
fprintf(' U   coarse      fine        texture     MCfull      +/-        EDoF        BCM(unif m)\n');
R=[];
for U = cfg.Uvec
    Pc=fama_floor('coarse',U,'hom',cfg.tau,cfg);
    Pf=fama_floor('fine'  ,U,'hom',cfg.tau,cfg);
    Pt=fama_floor('texture',U,'hom',cfg.tau,cfg);
    if quick, nt=2e4; else, nt=2e5; if U<=3, nt=1.2e6; end, end
    [Pm,ci]=fama_mc_floor('full',U,'hom',cfg.tau,nt,cfg);
    Pe=fama_bench('edof',U,cfg.tau,cfg); Pb=fama_bench('bcm',U,cfg.tau,cfg);
    fprintf(' %d  %.3e  %.3e  %.3e  %.3e  %.1e  %.3e  %.3e\n',U,Pc,Pf,Pt,Pm,ci,Pe,Pb);
    R=[R; U Pc Pf Pt Pm ci Pe Pb];
end

fprintf('\n--- Table III: anatomy of the gap ---\n');
fprintf(' U   m_I,b    c_b     marginal  coupling  total   measured A   Thm4/MC\n');
for k=1:size(R,1)
    U=R(k,1); mU=(U-1)*cfg.mI; cb=(mU+1)/(cfg.mI+1);
    [Ib,mIb]=fama_profile(U,'hom',cfg); Pm1=1;
    for b=1:cfg.B
        xb=cfg.mb(b)*Ib(b)*cfg.tau/(cfg.mb(b)*Ib(b)*cfg.tau+cb*cfg.Omb(b));
        Pm1=Pm1*betainc(xb,cfg.mb(b),cb);
    end
    fprintf(' %d  %7.2f  %6.3f  %8.2f  %8.2f  %6.2f  %10.2f  %7.2f\n', ...
        U,mIb(1),cb,Pm1/R(k,2),R(k,4)/Pm1,R(k,4)/R(k,2),R(k,5)/R(k,2),R(k,4)/R(k,5));
end

fprintf('\n--- Sec. IV-B: intra-block look count, U=6 ---\n');
for L=[1 2 4 8 16 25]
    c2=cfg; c2.Lb=L*ones(1,cfg.B);
    fprintf('  L_b=%2d : %.4e\n',L,fama_floor('fine',6,'hom',cfg.tau,c2));
end
fprintf('  Corollary 2 : %.4e\n',fama_floor('coarse',6,'hom',cfg.tau,cfg));

fprintf('\n--- Fig. 5: diversity slope ---\n');
ts=logspace(-1,-6,26);
for U=[3 5 8]
    P=arrayfun(@(t) fama_floor('coarse',U,'hom',t,cfg), ts);
    p=polyfit(log10(ts),log10(P),1);
    fprintf('  U=%d slope=%.4f   (d_inf=%.2f)\n',U,p(1),cfg.dinf);
end

fprintf('\n--- Fig. 6 and Sec. VII-G: hardening ---\n');
fprintf(' U  Ueff_hom Ueff_het  floor_hom   floor_het  ratio\n');
for U=cfg.Uvec
    [~,~,Uh]=fama_profile(U,'hom',cfg); [~,~,Ut]=fama_profile(U,'het',cfg);
    Ph=fama_floor('coarse',U,'hom',cfg.tau,cfg);
    Pt=fama_floor('coarse',U,'het',cfg.tau,cfg);
    fprintf(' %d  %7.2f %8.2f  %.4e  %.4e %6.2f\n',U,Uh,Ut,Ph,Pt,Pt/Ph);
end
fprintf('\n threshold sweep, het/hom ratio (sign reversal):\n tau(dB)   U=4      U=6      U=8\n');
for td=[-15 -10 -6 -3 0 3 6 10]
    t=10^(td/10); r=zeros(1,3); k=0;
    for U=[4 6 8], k=k+1; r(k)=fama_floor('coarse',U,'het',t,cfg)/fama_floor('coarse',U,'hom',t,cfg); end
    fprintf(' %6d  %7.3f  %7.3f  %7.3f\n',td,r(1),r(2),r(3));
end

fprintf('\n--- Table IV: selector comparison (Sec. VI) ---\n');
fprintf(' profile                 U  aware      agnostic   oracle     penalty\n');
P={[1 1 1 1],[1.6 1.2 0.7 0.5],[2.2 1.3 0.35 0.15]};
Nm={'flat [1 1 1 1]        ','moderate [1.6 1.2 .7 .5]','strong [2.2 1.3 .35 .15]'};
ntr = 2e4*(quick) + 1.5e5*(~quick);
for j=1:3
  for U=[4 6]
    o=fama_selector(U,P{j},25,ntr,cfg);
    fprintf(' %s %d  %.3e  %.3e  %.3e  %6.2f\n',Nm{j},U,o.aware,o.agnostic,o.oracle,o.agnostic/o.aware);
  end
end
save('results.mat','R','cfg');
end
