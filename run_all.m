%RUN_ALL  One command reproduces everything in the paper.
%
%   >> cd <this folder>
%   >> run_all
%
%   Produces
%     fama_log.txt          every number quoted in the paper (SEND THIS BACK)
%     fig2_sircdf.eps       Fig. 2
%     fig3_outage_snr.eps   Fig. 3   (slowest step, ~10-30 min)
%     fig4_floor_vs_U.eps   Fig. 4
%     fig5_diversity.eps    Fig. 5
%     fig6_hardening.eps    Fig. 6
%     fig3_data.mat, fig4_data.mat
%
%   Fig. 1 is a TikZ drawing inside main.tex and needs no MATLAB.
%
%   Set QUICK = true for a fast smoke test (wrong numbers, right shapes).

clear; close all; addpath(pwd);
QUICK = false;

diary off; if exist('fama_log.txt','file'), delete('fama_log.txt'); end
diary('fama_log.txt'); diary on;
fprintf('FAMA paper reproduction log   %s\n', datestr(now));
fprintf('MATLAB/Octave: %s\n\n', version);

%% ---- step 1: numerical self-test. Must pass before anything else. --------
fprintf('########## STEP 1: SELF-TEST ##########\n');
fama_selftest();

%% ---- step 2: every number quoted in Section VII -------------------------
fprintf('\n########## STEP 2: NUMBERS FOR SECTION VII ##########\n');
collect_results_core(QUICK);

%% ---- step 3: figures ----------------------------------------------------
fprintf('\n########## STEP 3: FIGURES ##########\n');
t=tic; fig2_sircdf();      fprintf('  Fig.2 done, %.0f s\n',toc(t));
t=tic; fig5_diversity();   fprintf('  Fig.5 done, %.0f s\n',toc(t));
t=tic; fig6_hardening();   fprintf('  Fig.6 done, %.0f s\n',toc(t));
t=tic; fig4_floor_vs_U();  fprintf('  Fig.4 done, %.0f s\n',toc(t));
t=tic; fig3_outage_snr();  fprintf('  Fig.3 done, %.0f s\n',toc(t));

fprintf('\n########## COMPLETE ##########\n');
fprintf('Send fama_log.txt back. The .eps files go next to main.tex.\n');
diary off;
