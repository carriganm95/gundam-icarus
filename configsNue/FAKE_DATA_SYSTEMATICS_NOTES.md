# Why do only the Template parameters pull in the fake-data studies?

Notes from investigating why, in the `nc`/`Q2`/`res0p5`/`res1p5` fake-data
fits, the `true_leading_electron_energy` Template parameters show large
post-fit pulls while the GENIE, Flux, and Detector systematics stay close to
their priors.

## TL;DR

Nothing is misconfigured or silently disabled. The behavior is explained by
two independent effects that compound each other:

1. **Structural**: the Template parameter set has no external covariance
   matrix, so GUNDAM automatically gives it a flat prior with **zero**
   penalty in the likelihood. It can move at no cost. The GENIE/Flux/Detector
   systematics *do* have a real Gaussian penalty (from their covariance
   matrices), so moving them costs chi2.
2. **Statistical**: the fake-data configs for `nc`/`Q2`/`res0p5`/`res1p5`
   normalize to ~2.36e19 POT, about 10.6x below the actual analysis target
   of 2.5e20 POT (used correctly in the `asimov` config). Lower POT weakens
   the statistical pull on *any* Gaussian-constrained parameter, making the
   systematics' already-small pulls even smaller.

Even at the correct POT, the Template will still dominate for anything that
falls inside its own subspace (see "Degrees of freedom" below) — that part
is expected/desired behavior for a template-based cross-section extraction.

## 1. What was actually checked

Rather than trust the config files alone, the fit *output* ROOT files
(`configsNue/data/<tag>/fit_<tag>.root`) were inspected directly
(`FitterEngine/preFit/parState_TNamed` vs `FitterEngine/postFit/parState_TNamed`)
to compare pre-fit vs. post-fit parameter values for each tag.

Result (values in normalized/sigma units, `useNormalizedFitSpace: true`):

| Parameter set              | # params | # that moved | max pull (Q2 tag) |
|-----------------------------|---------:|--------------:|-------------------:|
| Template (TrueEleE)         | 10       | 10             | ~0.9 sigma (using its display stdDev) |
| Cross section (GENIE)       | 52       | 43             | ~0.005 sigma |
| Flux (Beamline)             | 21       | 21             | ~0.005 sigma |
| Detector                    | 12       | 12             | ~0.013 sigma |

So the systematics genuinely *are* being fit — none are frozen/disabled —
they just land far closer to prior than the Template.

Of the 52 GENIE parameters, 9 (`ZExpPCAB1-4CCQE`, `MaNCEL`, `EtaNCEL`,
`DecayAngMEC`, `Theta_Delta2Npi`, `ThetaDelta2NRad`) show *exactly* zero
movement in every tag. Their dial branches were confirmed to exist correctly
in the input tree (`gundam_v27p2.root`), ruling out a naming mismatch (unlike
an old, stale `log.log` from March 2026 which *did* have a real
`ZExpPCAB1` vs. `ZExpPCAB1CCQE` naming bug — since fixed). The current
zero-movement pattern is consistent with `enablePca: true`
(`pcaDeltaChi2Threshold: 1E-4`) correctly freezing directions with
negligible chi2 response for this selection — worth a spot-check with
`enablePca: false` if it's ever in doubt, but not the cause of the
overall Template-vs-systematics asymmetry.

## 2. Root cause #1 (structural): Template has zero penalty cost

Confirmed directly in the GUNDAM source
(`GUNDAM-2.1.X/Repositories/gundam/src/ParametersManager/src/ParameterSet.cpp:1074-1084`):

```cpp
if( _priorFullCovarianceMatrix_ != nullptr ){
  par.setStdDevValue(std::sqrt((*_priorFullCovarianceMatrix_)[i][i]));
  par.setStepSize(std::sqrt((*_priorFullCovarianceMatrix_)[i][i]));
}
else{
  // stdDev will only be used for display purpose
  par.setStdDevValue(_nominalStepSize_);
  par.setStepSize(_nominalStepSize_);
  par.setPriorType(Parameter::PriorType::Flat);
  par.setIsFree(true);
}
```

`parameterSet_TrueEleE.yaml` has no `covarianceMatrixFilePath`, so every
Template parameter is auto-assigned `PriorType::Flat` and `isFree = true`.
The `nominalStepSize: 0.1` you set is explicitly a **display-only** value —
it is *not* used as a real prior width.

This is enforced in the likelihood itself
(`src/StatisticalInference/Likelihood/src/LikelihoodInterface.cpp:180-208`,
`evalPenaltyLikelihood`):

```cpp
double LikelihoodInterface::evalPenaltyLikelihood(const ParameterSet& parSet_) const{
  if( not parSet_.isEnabled() ){ return 0; }
  double buffer = 0;
  if( parSet_.getPriorCovarianceMatrix() != nullptr ){
      // ... buffer = delta^T * InverseCovarianceMatrix * delta ...
  }
  return buffer;
}
```

The Gaussian penalty term is only ever added `if getPriorCovarianceMatrix()
!= nullptr`. For the Template set that pointer is null, so `buffer = 0`
**unconditionally** (aside from the hard `[0,3]` bounds) — there is no
`1/sigma^2` term at all for Template, not even a "unit sigma" fallback.

Meanwhile GENIE/Flux/Detector each load a real covariance matrix
(`gundaminput_geniesyst.root`, etc.), giving `priorType = Gaussian`,
`stdDev = 1` (normalized fit space), confirmed by directly reading the
`priorType_TNamed`/`stdDev_TNamed` entries in the propagator dump of
`fit_Q2.root`.

This is the standard, correct design for a template/unfolding-style
cross-section extraction: the Template parameters *are* the measurement
(cross section per true bin) and should have no external Bayesian prior, so
as not to bias the extraction. The systematics represent uncertain inputs
with real external constraints (beam simulation, hadron production data,
calibration), which is exactly what justifies their Gaussian penalty.

## 3. Root cause #2 (statistical): why POT scale changes the pull size

### What `chi2_stat` actually is in GUNDAM

Your fitter config selects `jointProbabilityConfig: type: "PoissonLLH"`.
This picks GUNDAM's `PoissonLogLikelihood::eval`
(`src/StatisticalInference/JointProbability/include/PoissonLogLikelihood.h`),
which computes, **per bin**, the saturated Poisson deviance:

```cpp
llhValue = 2.0 * (predVal - dataVal + dataVal * TMath::Log(dataVal / predVal));
```

i.e. for bin `i` with MC prediction `mu_i` and observed/fake-data count `n_i`:

```
D_i = 2 * ( mu_i - n_i + n_i * ln(n_i / mu_i) )
```

This is the standard "-2 * ln(likelihood ratio)" statistic: it compares the
Poisson likelihood of observing `n_i` given prediction `mu_i` against the
best-possible ("saturated") model where the prediction exactly equals the
data. It is exactly zero when `mu_i = n_i`, and grows the further `mu_i`
strays from `n_i`.

Two special cases handled explicitly in the code:
- `mu_i <= 0` (a bin with zero MC events but nonzero data expectation)
  returns `+infinity` for that bin. This is exactly what produced the
  `Reference stat log-likelihood for PCA: inf` line seen in the old, stale
  `log.log` from March 2026 — a sign that some bin had zero simulated MC
  weight while being compared against a nonzero expectation.
- `n_i <= 0` uses the `x*ln(x) -> 0` limit, reducing to `llh_i = 2*mu_i`.

`chi2_stat` (what earlier sections in this note call `chi2_stat(theta)`) is
just the sum of `D_i` over every bin, in every sample, in every dataset:

```cpp
// LikelihoodInterface::evalStatLikelihood
statLikelihood = sum over samples of jointProbability->eval(samplePair)
// jointProbability->eval sums D_i over every bin in that sample
```

and the total objective GUNDAM actually minimizes is:

```cpp
// LikelihoodInterface::evalLikelihood
totalLikelihood = statLikelihood + penaltyLikelihood
```

— i.e. exactly `chi2_stat(theta) + penalty(theta)`, the same structure used
throughout this note. `penaltyLikelihood` is the sum of the
`evalPenaltyLikelihood(parSet)` terms from Section 2 (zero for Template,
`delta^T * InvCov * delta` for GENIE/Flux/Detector).

### Why `D_i` behaves like a chi2 (and why "Fisher information" is the right
language)

In the high-statistics / small-mismatch limit (`mu_i` close to `n_i`), a
Taylor expansion of `D_i` recovers the familiar chi2 form. Writing
`mu_i = n_i + delta` for small `delta`:

```
D_i = 2*( delta + n_i * ln(n_i/(n_i+delta)) )
    ~= 2*( delta - delta + delta^2/(2*n_i) )
    =  delta^2 / n_i
    =  (mu_i - n_i)^2 / n_i
```

which is exactly the standard Pearson chi2 term. This confirms that the
Poisson deviance GUNDAM actually computes reduces to a normal chi2 once
counts are reasonably large, and that its curvature with respect to a
parameter `theta` (via `mu_i(theta)`) is the Fisher information
`I ~ sum_i (dmu_i/dtheta)^2 / n_i` used below — scaling linearly with the
overall count/POT level, since `mu_i` and `n_i` both scale linearly with
`c` while `dmu_i/dtheta` scales with `mu_i` (hence with `c`) too.

### The likelihood balance

Each Gaussian-constrained parameter `theta` (prior 0, width `sigma`) is fit
by minimizing:

```
chi2(theta) = chi2_stat(theta) + (theta - 0)^2 / sigma^2
```

Near the minimum, `chi2_stat(theta)` is well approximated by a parabola:

```
chi2_stat(theta) ~= I * (theta - theta_hat)^2
```

where `theta_hat` is what the data alone would prefer (ignoring the prior),
and `I` is the Fisher information — the curvature of that parabola, i.e.
how sharply the data constrains `theta`.

### I scales linearly with POT

For a Poisson process, `I ~ sum (dmu/dtheta)^2 / mu`. Scaling every
predicted/observed count by a constant POT factor `c` (`mu -> c*mu`) gives:

```
I = c * I_0
```

i.e. **Fisher information scales linearly with your POT/statistics scale.**
More data -> more statistical leverage on every parameter, systematics and
Template alike.

### Solving for the equilibrium pull

Minimizing the total chi2 gives:

```
theta* = theta_hat * r / (1 + r),        r = I * sigma^2
```

`r` is the ratio of "how much the data knows about theta" to "how much the
prior already constrains it":

- `r << 1` (low POT, or a tight/informative prior): `theta* ~= theta_hat * r`
  — the pull is suppressed, and scales *linearly* with POT in this regime.
- `r >> 1` (high POT, or a very loose/absent prior): `theta* -> theta_hat`
  — the data fully wins, prior is irrelevant.

This is why cutting POT ~10.6x (2.5e20 -> 2.36e19) shrank the systematic
pulls by roughly the same order — not a bug, just `r` scaling down linearly
with `c`.

### Spring analogy

Think of two springs pulling on `theta`:
- A "data spring" anchored at `theta_hat`, stiffness `I` (grows with POT).
- A "prior spring" anchored at 0, fixed stiffness `1/sigma^2` (set by the
  external covariance matrix, independent of your dataset).

At low POT the data spring is floppy and the stiff prior spring wins —
`theta` stays near 0. At high POT the data spring stiffens until it's worth
paying the prior penalty to relieve the (now much more demanding)
statistical mismatch, so `theta` moves toward `theta_hat`.

### Why this doesn't apply to the Template the same way

For Template, `sigma^2 -> infinity` (no covariance matrix => no prior
spring at all, confirmed in Section 2). Plug that into `r = I * sigma^2`:
`r -> infinity` for *any* nonzero `I`, so `theta* -> theta_hat` regardless
of POT. There's no opposing spring to contest, so it doesn't need to "win" a
stiffness contest — it moves fully to `theta_hat` even with very little
statistical leverage behind it.

What POT *does* still control for the Template is not *whether* it moves,
but how precisely `theta_hat` itself is determined (i.e. the post-fit
uncertainty on the extracted value) — lower POT means more Poisson-noise
smearing around the same fully-migrated central value, not a smaller pull.

## 4. Why doesn't the Template just absorb the *entire* discrepancy?

Because "free to move" only means it takes 100% of whatever it is *capable*
of explaining — not that it has the degrees of freedom to explain
everything.

The Template is one number per true-electron-energy bin
(`binning_TrueEleE.txt`, ~10 bins), applied uniformly to every event
originating in that true bin. That is the *entire* subspace of predictions
it can reach. Any part of a data/MC discrepancy that is expressible purely
as "rescale the rate from true-E bin i" gets absorbed completely, for free.
Anything outside that subspace is simply unreachable by Template, no matter
how unconstrained it is:

- The **`nc`** fake-data formula reweights by
  `true_category == 4 ? f(reco_leading_electron_energy) : 1.0` — a function
  of **reco** electron energy. Template can rescale the *total* rate from a
  true bin, but cannot reshape how those events are distributed across reco
  bins after smearing/reconstruction. That reco-level shape component is
  invisible to Template and is exactly the kind of effect a
  detector-energy-scale/smearing systematic could, in principle, reach.
- The **`Q2`** formula depends on `true_momentum_transfer` and `category`,
  a different true-level variable than the Template's binning axis. Only
  the part of that effect that happens to project onto true electron energy
  gets absorbed "for free"; the rest is outside Template's reach.
- The cosmic (`...Cosmics`) and beam-off backgrounds don't have a physical
  "true leading electron energy" at all, so Template (built on that truth
  variable) structurally cannot correct anything in those background
  components.

So: **Template affects reco bins only through a uniform rescaling of a true
bin's total contribution, propagated through the (fixed) smearing/migration
matrix — it cannot independently reshape the reco-level distribution from a
given true bin.** The residual piece of any injected mismodeling that lies
outside the "true-E-bin-normalization" subspace is what's actually available
for the Gaussian-constrained systematics (or, if no dial spans that
direction, for irreducible residual chi2) to respond to. This is the origin
of the small but nonzero GENIE/Flux/Detector pulls observed (~1e-3 to 2e-2
sigma) — not statistical noise, but the genuinely non-degenerate component
of the injected effect.

## 5. What this means for the stated goal of the fake-data studies

Stated goal (from the analysis write-up): verify (1) GUNDAM can recover a
cross section that differs from the nominal GENIE AR23_20i prediction, and
(2) that residual mismodeling not covered by a direct nuisance knob does not
significantly bias the extracted result.

Given the mechanisms above:

- **Goal 1** is tested by comparing the **post-fit Template values against
  the true injected per-true-bin scaling** (known exactly, since the fake
  data is built from a specific reweight formula). This is the decisive
  check — not whether the systematics moved.
- **Goal 2** is tested by looking for a **systematic offset between the
  recovered Template values and the true injected values**, and/or
  significant residual post-fit chi2/tension. Small systematic pulls near
  prior are not, by themselves, evidence of a problem — they can equally
  mean the injected mismodeling was fully and correctly captured by the
  Template (the intended behavior), with only a small, correctly-identified
  non-degenerate residual left for the systematics.

A **different** question — "how much of this mismodeling could my
systematic model explain if the Template weren't available to absorb it
first" — is a distinct diagnostic (systematic-model sufficiency/headroom),
not what the two stated goals ask for. It would require re-running with the
Template parameter set disabled (`isEnabled: false` in
`parameterSet_TrueEleE.yaml`) so only GENIE/Flux/Detector can respond, then
comparing pulls and residual chi2. Useful as a separate check, but not
necessary to answer goals 1 and 2 above.

## 6. Action items

- [ ] Align `nc`/`Q2`/`res0p5`/`res1p5` fake-data POT normalization with the
      actual analysis exposure (2.5e20, matching `asimov`), instead of the
      current ~2.36e19, so pull sizes reflect real expected sensitivity
      rather than an artificially statistics-starved test.
- [ ] For each tag, compute and compare post-fit Template values against the
      true injected per-true-bin scaling factor implied by the fake-data
      reweight formula, to directly quantify closure (goal 1) and bias
      (goal 2).
- [ ] (Optional, separate diagnostic) Re-run with the Template parameter set
      disabled to gauge how much of each injected effect the systematic
      model could explain on its own.
