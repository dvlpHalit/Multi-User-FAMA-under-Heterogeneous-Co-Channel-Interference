# Outage Analysis of Multi-User FAMA under Heterogeneous Interference over Block-Correlated Channels with Nakagami-m Textures
Reproduction code for the outage-floor analysis of a fluid antenna terminal receiving
heterogeneous co-channel interference over block-correlated Nakagami-m channels.
> \*\*Paper.\*\* M. H. Aldemir, M. A. Karabulut, and A. F. M. Shahen Shah, "Multi-User FAMA
> under Heterogeneous Co-Channel Interference: Outage Floors over Block-Correlated
> Nakagami-\*m\* Channels," submitted to \*IEEE Transactions on Communications\*, 2026.
Every figure and table in the paper is produced by a script in this repository. Nothing is
hand-tuned and there are no fitted constants anywhere in the analytical path.
---
Scope
The code covers four things:
The interference reduction. The second-order (Welch–Satterthwaite) collapse of a
heterogeneous multi-source interference sum to a single Gamma law, and the resulting
Fisher–Snedecor block SIR (Lemmas 1 and 2).
The closed-form floors. The fine floor of Theorem 1, the single-look floor of
Corollary 2, the hardened floor of Corollary 3, and the diversity order of Theorem 2.
The texture-aware correction. The parameter-free calibration of Proposition 2 and
the corrected floor of Theorem 4, together with the FGM structural ceiling of Theorem 3.
The correlated-field simulation. Monte Carlo over the physical Abdi–Kaveh
correlation, which is what the analytical expressions are validated against.
The analysis is validated against the physical correlation model, not against the banded
surrogate it approximates. This is deliberate: the banded law in Eq. (2) is a structural
device used to expose the geometry behind the effective rank and the block-boundary
discussion, and is never substituted numerically for the physical kernel.
---
Requirements
	
MATLAB	R2021a or later
GNU Octave	8.4 or later
Toolboxes	none required
Every analytical routine uses base functions only (`betainc`, `gammainc`, `gamma`,
`besseli`, `chol`). This is not incidental. The closed form of Theorem 1 contains Tricomi's
confluent hypergeometric function, which in MATLAB lives in `kummerU` (Symbolic Math
Toolbox) and has no base-Octave equivalent. The numerical path deliberately avoids it by
evaluating the pre-expansion form of Remark 4 instead, so the repository runs on a stock
installation. Eq. (16) is kept in the paper for interpretation, not for evaluation.
Gamma variates are drawn with a local Marsaglia–Tsang sampler so that results are identical
across MATLAB and Octave for a given run index.
---
Quick start
```matlab
cd src
run\_all                 % regenerates every figure and table, \~N minutes
```
or, for a single result,
```matlab
fig4\_floor\_vs\_U         % writes fig4\_floor\_vs\_U.eps and prints Table IV
```
Figures are written as `.eps` at the sizes used in the manuscript. Set
`cfg.mc\_trials = 1e4` in `config.m` for a fast smoke test; the published numbers use
`1.2e6`.
---
Layout
> Rename to match your working files if these differ.
```
src/
  config.m                  % Table III: the single source of truth for all parameters
  build\_J.m                 % Abdi-Kaveh von Mises correlation, Eq. (1)
  effective\_rank.m          % Eq. (3); fixes B and L\_b from the aperture
  interference\_reduction.m  % Lemma 1 and Definition 1: Ibar\_b, m\_I,b, U\_eff,b
  block\_sir\_cdf.m           % Lemma 2: scaled Fisher-Snedecor CDF, Eq. (11)
  floor\_theorem1.m          % Eq. (15) via the cancellation-free route of Remark 4
  floor\_corollary2.m        % Eq. (17), single-look
  floor\_corollary3.m        % Eq. (19), hardened
  floor\_theorem4.m          % Eq. (25), texture-aware
  texture\_calibration.m     % Proposition 2: m\_upsilon and c\_b, no free parameters
  fgm\_ceiling.m             % Theorem 3: the 2^(B-1) bound
  mc\_correlated\_field.m     % Monte Carlo over the physical field
  selectors.m               % aware / agnostic / oracle rules of Table II
  run\_all.m

  fig2\_sircdf.m             % Fig. 2
  fig3\_outage\_snr.m         % Fig. 3
  fig4\_floor\_vs\_U.m         % Fig. 4(a,b) and Table IV
  fig5\_diversity.m          % Fig. 5
  fig6\_hardening.m          % Fig. 6
  table2\_selectors.m        % Table II

ref/
  reference\_values.json     % high-precision values for regression testing
  verify.m                  % checks every script against ref/ to 1e-9
```
---
What produces what
Paper object	Script	Kind
Table III (configuration)	`config.m`	parameters
Fig. 2, block SIR CDF	`fig2\_sircdf.m`	analysis + MC
Fig. 3, outage vs pilot SNR	`fig3\_outage\_snr.m`	MC + Theorem 4
Fig. 4(a), floor vs U	`fig4\_floor\_vs\_U.m`	Thms 1, 4; Cor. 2; MC
Fig. 4(b), required amplification vs FGM ceiling	`fig4\_floor\_vs\_U.m`	analysis
Fig. 5, floor vs threshold	`fig5\_diversity.m`	analysis only
Fig. 6, hardening at equal mean power	`fig6\_hardening.m`	analysis + MC
Table II, three selectors	`table2\_selectors.m`	MC
Table IV, anatomy of the residual	`fig4\_floor\_vs\_U.m`	analysis + MC
Fig. 5 is analytical throughout. Floors below `1e-10` are beyond Monte Carlo reach, so that
figure verifies the exponent of the model rather than of the correlated field.
---
Reproduction targets
Analytical values are deterministic and should match to the digits shown. Monte Carlo values
carry sampling error; the deepest points are the least constrained.
Configuration (`config.m`): `N = 100`, `W = 3` wavelengths, `dx = 0.030` wavelengths,
`kappa\_A = 8`, `B = 4`, `L\_b = 25`, `|S\_b| = 8`, `m\_b = \[3.0 2.4 1.8 1.4]`, `m\_I = 2.5`,
per-emitter SIR `10 dB`, threshold `tau = 0 dB`, path-loss exponent `3.5`.
Structural quantities
Quantity	Value
Aperture effective rank, Eq. (3) on Eq. (1)	`4.07`
Per-block effective rank, `L\_b = 25`	`1.399`
Decay ratio matching the per-block rank	`rho\_b = 0.945`
PSD limit of the banded law at `d = 4`	`rho\_b <= 0.6946`
Floor diversity order	`d\_inf = 8.6`
Floors at `tau = 0 dB` (analytical, `tau\_S = 1`)
U	`m\_I,b`	`c\_b`	Corollary 2	Theorem 4	MC field floor
2	2.5	1.000	`1.412e-07`	`1.886e-05`	`7.5e-06`
3	5.0	1.714	`1.403e-05`	`1.951e-04`	`1.36e-04`
4	7.5	2.429	`1.962e-04`	`9.494e-04`	`5.40e-04`
5	10.0	3.143	`1.163e-03`	`3.092e-03`	`1.83e-03`
6	12.5	3.857	`4.255e-03`	`7.781e-03`	`4.09e-03`
8	17.5	5.286	`2.480e-02`	`3.001e-02`	`1.36e-02`
The ratio Theorem 4 / MC is `2.52, 1.44, 1.76, 1.69, 1.91, 2.20`. It is one-signed at every
emitter count, so it is a bias and not sampling noise, and it lies in the conservative
direction.
Consistency checks that must pass
(S1) `Ibar\_b -> 0` drives Eq. (15) to zero identically, to machine precision.
(S2) With homogeneous emitters the beta argument `x\_b = m\_b/(m\_b + 25)` is invariant to
U; the whole dependence on U enters through the second shape parameter.
(S3) `kappa\_A -> 0` collapses Eq. (1) to the Clarke kernel `besselj(0, D)`.
The equivalent-degrees-of-freedom Rayleigh floor coincides with Theorem 1 at `B = 1`,
`L\_1 = 2\*ceil(W)+1`, `m\_b -> Inf` to thirteen significant digits.
The diversity slope, fitted below `tau\_S = 1e-4`, returns `8.6` to five decimals.
---
Numerical notes
Three things will silently produce wrong answers if changed.
1. The alternating sum in Eq. (15) must not be evaluated directly. For the configuration
above, the largest summand exceeds the block value by twelve orders of magnitude at
`tau\_S = 1e-1` and by eighteen at `tau\_S = 1e-3`. Double precision is exhausted well before
the deep-floor regime. `floor\_theorem1.m` evaluates the pre-expansion form of Remark 4
instead, where both integrands lie in `\[0,1]` and no large quantities are subtracted. The
two routes agree to `4e-10` where the closed form is still evaluable, and the
cancellation-free route agrees with high-precision references to `1e-15` throughout.
2. Gauss–Laguerre quadrature is unsuitable for the texture integral. The kernel
`(1 - exp(-a/t))^L` concentrates its mass near `t -> 0+`, where Gauss–Laguerre nodes have no
support. Using it for the inner expectation underestimates the block factor by roughly two
orders of magnitude at `m\_b = 1`, `L\_b = 16`. The inner expectation is therefore taken by
the inverse-CDF substitution `t = F\_t^{-1}(p)` on a panel grid anchored at the half-power
texture; generalized Gauss–Laguerre is used only for the outer expectation over the
interference, where the integrand is smooth.
3. The correlation law is the exponent, not a scale. In `build\_J.m` the angular
parameters enter inside the Bessel argument. Applying any severity-dependent factor as an
external multiplier instead makes the effective rank independent of m, which quietly
breaks the partition rule and every floor level that follows from it.
Two further points worth knowing. The block covariance is mathematically full rank but its
eigenvalues fall off as `20.76, 3.92, 0.30, 1.2e-2, 3e-4` and then to machine zero, so the
diversity asymptote is entered late; this is why the slope fit tightens only below
`tau\_S = 1e-4`. And the physical correlation at one block spacing is still `0.241`, with
cross-block port pairs carrying `27%` of the total correlation energy, which is the
mechanism behind the residual that Theorem 4 repairs.
---
Citation
```bibtex
@article{Aldemir2026FAMAfloor,
  author  = {Aldemir, Mevlana Halit and Karabulut, Muhammet Ali and
             Shahen Shah, A. F. M.},
  title   = {Multi-User {FAMA} under Heterogeneous Co-Channel Interference:
             Outage Floors over Block-Correlated {Nakagami-$m$} Channels},
  journal = {IEEE Transactions on Communications},
  year    = {2026},
  note    = {Submitted}
}
```
Please update this entry once the paper is accepted.
---
Funding
Supported by the Scientific and Technological Research Council of Türkiye (TÜBİTAK) under
Project 124E519.
License
Released under the MIT License. See `LICENSE`.
Contact
Mevlana Halit Aldemir — `mevlanahalit.a@msu.edu.tr`
Department of Electronic Communication Technologies, National Defence University, Balıkesir,
Türkiye.
Issues and pull requests are welcome. If a number in the reproduction table does not come out
as listed, please open an issue with your MATLAB or Octave version and the output of
`ref/verify.m`.
