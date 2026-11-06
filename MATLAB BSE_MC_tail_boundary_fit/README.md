# BSE and MC chain closures

This package compares five closures of the chain representation of a structured
bosonic environment. The chain is kept exactly up to site N, and the closure
replaces everything beyond site N:

1. **Asymptotic BSE:** `kappa(N)^2 / (p + 1i*omega(N+1) + SigmaInf(p))`.
2. **Chain extended BSE:** the same form with `SigmaBar_N(p)`, whose parameters are
   obtained from the available chain coefficients.
3. **Semicircle-fit MC:** the published Markovian closure with M damped modes.
4. **MC: chi fit to SigmaBar:** an MC-form network with M damped modes, fitted by the
   cut-invariant error to the chain-extended BSE closure at site N.
5. **MC: chi fit to exact tail:** the same network, fitted to the exact tail
   `Sigma_exact_(N+1 -> N)`. This closure needs J(omega).

Every MC closure keeps the exact chain up to site N and couples its M damped modes
to site N, so it has N chain modes plus M damped modes. The default is M=10 and
N=10,...,60.

## Run

Open this folder in MATLAB and run:

```matlab
results = run_structured_bse_benchmark;
```

The package includes the universal chi fit of the M=10 network to the semicircle,
which the default run reuses when its settings match. The run fits the two
chi-fitted MC closures at every N, computes the reference and closure dynamics,
the errors, the tables, and nine figures. The chi fits at the 51 values of N take
most of the run time. All outputs go to `results/`. Base MATLAB is sufficient; no
toolbox is needed for any fit.

```matlab
run_structured_bse_benchmark('plot');   % regenerate the nine figures
report = verify_asymptotic_benchmark;   % independent numerical checks
```

`'fit'` computes and saves the chain coefficients and the MC closures without
dynamics; `'evaluate'` uses saved fits with matching settings. Edit
`benchmark_settings.m`, or pass overrides as the second argument, for example:

```matlab
run_structured_bse_benchmark('all', struct('mcPortLink', 'kappaInf'));
```

## Closures and settings

- **Asymptotic BSE:** `omegaInf=0.5` and `kappaInf=0.25`. Site N+1 is eliminated
  analytically with its own coefficient `omega(N+1)` and coupling `kappa(N)`.
- **Chain extended BSE:** the same form with `omegaBar_N=omegaInf+a/(N+1)^2` and
  `kappaBar_N=kappaInf+b/(N+1)^2`, where a and b are fitted to the last third of
  `omega(1:N+1)` and `kappa(1:N)` (`chain_only_bse_parameters.m`). No exact tail is
  fitted.
- **Semicircle-fit MC:** the full-precision published universal coefficients for
  M=10 (`mc_coefficients/`), fitted to the correlation function of the Wigner
  semicircle. Internal frequencies, hoppings and losses are shifted and rescaled with
  `omegaInf` and `kappaInf`. The couplings to site N are `kappa(N)*c`
  (`mcPortLink='kappaN'`), so this closure approximates
  `kappa(N)^2 / (p + 1i*omegaInf + SigmaInf(p))`. `mcPortLink='kappaInf'` uses
  `kappaInf*c`, the rescaling rule of the published MC.
- **MC: chi fit to SigmaBar:** the M=10 network is first fitted once to the semicircle
  by minimizing the cut-invariant error, starting from the published coefficients
  (`universal_chi_mc_networks.m`, cached in `cache/chi_mc_universal.mat`). At each N,
  all network parameters are then fitted to
  `kappa(N)^2 / (p + 1i*omega(N+1) + SigmaBar_N(p))` on the band interior of
  `SigmaBar_N` (`chi_mc_model.m`). This closure uses only `omega(1:N+1)` and
  `kappa(1:N)`.
- **MC: chi fit to exact tail:** at each N, all network parameters are fitted to the
  exact tail in the band, starting from the previous closure at the same N.
- Both chi fits add a weak out-of-band penalty with the chain-extended BSE closure as
  its target, since the inverse recurrence for the exact tail is unstable outside the
  band. The fits use a base-MATLAB Levenberg-Marquardt routine (`lm_least_squares.m`)
  with an analytic Jacobian (`mc_network_response.m`).

The structured Ohmic bath uses `alpha=0.02`, Gaussian amplitude `4`, center `0.70`,
and width `0.07`, on `[0,1]`. The qubit frequency is `0.5`. The benchmark uses
`T=2000`, population time step `0.1`, and spectral broadening `etaSpec=1/T=0.0005`.
The continuum reference is evaluated directly from the spectral density, including
bound-state contributions when present. The BSE time-nonlocal response is not an
additional explicit oscillator; equal N is not a claim of equal computational cost.

## Error measures

- **Relative population error** (Figure 3):
  `epsP = sqrt( int_0^T |P_X - P_ref|^2 dt / int_0^T |P_ref|^2 dt )`.
  Writing `A_X = A + dA`, the error `P_X - P = 2 Re(conj(A) dA) + |dA|^2` is dominated by
  the cross term, which depends on the phase of the reflected amplitude dA relative to
  A. This phase advances with the round trip to site N, so epsP oscillates with N, with
  a period of about 12 sites here. The CSV files also list the maximum population error,
  the relative amplitude error `sqrt( int_0^T |dA|^2 dt / int_0^T |A|^2 dt )`, and the
  phase-independent bound
  `sqrt( int_0^T (2|A||dA| + |dA|^2)^2 dt / int_0^T |P_ref|^2 dt ) >= epsP`; both
  decrease smoothly with N.
- **Relative spectral error** (Figure 4): the same normalized L2 error of
  `real(Sigma_port(etaSpec - 1i*omega))` over the full band `0 <= omega <= 1`, the range
  over which the published MC is compared with the semicircle (Nusseler et al. 2022,
  Fig. 2(c)). This absolute error stays finite at the band edges, so it includes them.
- **Cut-invariant error** between two self-energies presented to the same site,

  ```
  chi(omega) = |Sigma_a - Sigma_b| / |Sigma_a + conj(Sigma_b)|,   p = 0+ - i*omega.
  ```

  On the real axis every step of the self-energy recurrence leaves chi unchanged, so
  chi at the cut equals chi between the two responses at the system port. chi=1 wherever
  one of the two self-energies has zero real part, which happens at the band edges for
  any closure, so `chi_rms` is its root mean square over the band interior
  `|omega-omegaInf| <= chiBandFraction*2*kappaInf` (default 0.98, so
  `0.01 <= omega <= 0.99`). The chi fits use the same window. Figure 9 shows chi over the
  full band.
  - **Total error** `chi_tot` (Figures 5 and 9): each closure against the exact tail
    `Sigma_exact_(N+1 -> N)`. For an MC closure it is bounded by its closure and depth
    errors, `|chi_M-chi_N|/(1-chi_M*chi_N) <= chi_tot <= (chi_M+chi_N)/(1+chi_M*chi_N)`,
    so once the depth error of its target falls well below its closure error, `chi_tot`
    is close to the closure error. This floor is set by the M damped modes and does not
    decrease with N.

  As in Sec. III, the total error is split into two parts:
  - **Depth error** `chi_N` (Figure 6): each BSE closure at site N against the exact tail
    `Sigma_exact_(N+1 -> N)`.
  - **Closure error** `chi_M` (Figure 7): each MC closure against its own target at site
    N. The semicircle-fit MC has the asymptotic BSE closure as its target, the chi fit to
    SigmaBar the chain-extended BSE closure, and the chi fit to the exact tail the exact
    tail, so its closure error is also its total error.
  For the BSE closures the total error equals the depth error. The totals of all
  closures are listed in `closure_errors_vs_N.csv` (`ChiRms`) and `chi_budget_vs_N.csv`.

The exact tail is obtained from the exact system-port response `J + i*Delta` by the
inverse recurrence (`prepare_mc_chi_targets.m`). Inside the band each inverse step
preserves chi, so this transport is stable.

## Figures

1. `Figure1_time_domain`: populations at N=15 on `[0,100]`, with the exact continuum
   and the hard cutoff.
2. `Figure2_effective_spectrum`: broadened port spectra at N=15, with zoomed panels of
   the peak and the upper band edge. Dotted boxes in the main panel mark the zoomed
   windows.
3. `Figure3_population_error_vs_N`: relative population error for N=10,...,60.
4. `Figure4_spectral_error_vs_N`: relative spectral error for N=10,...,60.
5. `Figure5_total_error_vs_N`: total error `chi_tot` (rms) of the five closures against
   the exact tail for N=10,...,60.
6. `Figure6_depth_error_vs_N`: depth error `chi_N` (rms) of the two BSE closures for
   N=10,...,60.
7. `Figure7_closure_error_vs_N`: closure error `chi_M` (rms) of the three MC closures for
   N=10,...,60.
8. `Figure8_self_energy_at_cut`: `Re Sigma_(N+1 -> N)(0+ - i*omega)` of the three MC
   closures and their targets at N=50 (`Ncut`), with an inset of the upper band edge.
9. `Figure9_total_error_at_cut`: total error `chi_tot(omega)` of the five closures against
   the exact tail at N=50. The exact tail is the reference and has no curve.

Legends are inside the axes; Figures 3, 4, 5, 8 and 9 list the BSE closures in the first
column and the MC closures in the second. Legend text uses font size 17, close to the
caption size when a figure is reduced to half its width in the manuscript. In Figures 3-9,
which are panels in the manuscript, each legend is centered and at most 75% of the axes
width, so the top-left corner stays free for the panel label; the font is reduced until
the legend fits. The upper y-limit of each figure
is raised until the curves under the legend lie below it. Each closure keeps one colour
and line style in all figures. Each figure is saved as PDF, PNG, and editable MATLAB FIG.

## Tables, data, and checks

- `results/Table_I.tex`: chain-extended BSE parameters and relative population errors at
  N=15, 30, 45, 60 (`Ntable`).
- `results/closure_errors_vs_N.csv`: Figures 3-5 data, maximum population errors,
  population error bounds, relative amplitude errors, evaluation times, and the
  chain-extension parameters at every N.
- `results/chi_budget_vs_N.csv`: Figures 6-7 data, the total `chi_rms` of the MC
  closures, and the error of the published MC against the response it approximates.
- `results/time_domain_curves.csv`, `spectral_curves.csv`: curves at N=15, 30, 45, 60,
  with the exact reference and the hard cutoff.
- `results/frequency_comparison_at_N50.csv`: Figures 8 and 9 data.
- `results/mc_closure_networks.csv`: physical parameters of the three MC closures at
  every N (Closure 1: semicircle fit, 2: chi fit to SigmaBar, 3: chi fit to exact tail).
- `results/chi_mc_universal_network.csv`: the universal chi fit in units of kappaInf
  (mode frequencies `omegaInf+kappaInf*h`, hoppings `kappaInf*t`, rates
  `kappaInf*gamma`).
- `results/benchmark_results.mat`, `fitted_environments.mat`: reproducible data.
- `results/verification_asymptotic.json`: closure stability and passivity, independent
  matrix-exponential propagation, refined time, frequency, and real-axis quadratures
  for every N and closure, refined chain coefficients, cut invariance of chi, the chi
  budget, the population error bound, and the Figure 3 data.

## Changes from the previous release

- CAP is removed from the code, the cache, and all figures.
- The MC closures no longer split N into chain modes and damped modes. They use N chain
  modes plus M=10 damped modes attached to site N, so no allocation is selected.
- The published MC couples to site N with `kappa(N)` by default.
- The targets at site N include `omega(N+1)` and `kappa(N)`, and the two chi-fitted MC
  closures are fitted at every N.
- The chi fits and `chi_rms` use the band interior `0.01 <= omega <= 0.99` (previously
  `0.02 <= omega <= 0.98`). The relative spectral error uses the full band
  `0 <= omega <= 1`, as in the comparison of the published MC with the semicircle.
- Figure 3 uses the relative population error, and Figures 3-7 scan N=10,...,60.
  The previous Figures 5 and 7 are merged into Figure 3. The cut-invariant error is
  shown as the total error of all five closures (Figure 5), the depth error of the BSE
  closures (Figure 6) and the closure error of the MC closures (Figure 7). Figure 8 shows
  the MC closures and their targets at N=50, and Figure 9 the total error at N=50.
