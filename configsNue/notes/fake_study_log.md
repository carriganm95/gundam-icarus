# Fake-study run log — ICARUS NuMI nueCC GUNDAM

One row per study run: what changed and the headline result. **Add a row for
every new run** (newest at the bottom of its section).

- "Mean xsec unc" = mean over bins of the fractional error on the unfolded
  cross section (`calcXsec/histograms/*_XS`), computed from `xsec_*.root`.
- Unless noted, runs use `run_fake_study.sh <tag>` with the default 2.5e20 POT
  and NTOYS=1000.
- Entries before 2026-10-27 were reconstructed from output folders, run logs and
  the config stored in each output ROOT file. Toy counts and POT for those runs
  weren't checked.

## Fake-data definitions (what each tag changes)

Each tag differs from `asimov` only in how the fake data is made (from
`fakeStudyConfigs/config_Fitter_FakeData_<var>_<tag>.yaml`). The MC model is
unchanged.

| Tag | Fake-data change |
|---|---|
| `asimov` | None (fake data = nominal MC) |
| `nc` | NC background (`true_category==4`) scaled ×1.3 / 1.15 / 0.95 / 0.85 for reco E_e <0.5 / <0.8 / <1.2 / ≥1.2 GeV |
| `Q2` | Weight suppressing low momentum transfer: up to −70% in 0 < q < 0.7, quadratic shape |
| `res0p5` / `res1p5` | RES events (`true_interaction_mode==1`) ×0.5 / ×1.5 |
| `pAngle` | Backgrounds (`true_category!=0`) weighted by reco p_e·sinθ_e (GeV) |
| `emb` | Fake data from the `recombination` detector-variation tree (EMB recombination), `merged_final_v31e_flux.root` |
| `gainHigh` / `gainLow` | Fake data from the `induction_high` / `induction_low` detector-variation trees |
| `data_10pct` | 10% real-data sample (`run_data_10pct.sh`), not a fake study |

Systematic-breakdown studies (from `run_frac_uncertainty.sh`, Asimov):
`frac_uncertainty_asimov` (systematic subsets: all / only / N−1),
`_diagonal` (`--diagonal`), `_noflux` (`--zero-flux-norm`),
`decompose_asimov` (`--decompose`: stat/template/MC-throw/flux-norm split),
`split_detector_asimov[_nomcthrow]` and `split_genie_asimov` (`--split`, one dial at a time).

## Run history

### Early runs (Feb–Sep 2026), single output files in `data/EleE/`

| Date | Output | Change / purpose |
|---|---|---|
| 2026-02-26 → 03-17 | `example*.root`, `2026_new_CV_fitter_output_wpca.root` | Workshop examples / first CV fit with PCA |
| 2026-04-17 → 04-20 | `nue_asimov*.root`, `nue.root`, `nue_XSec.root`, `nue_asimovXSec.root` | First nue Asimov fit + xsec, flat weight 0.321 (no PPFX) |
| 2026-05-15 | `nue_fitter_QE0p2/0p7/1p2.root` | Fake data: QE (`true_interaction_mode==0`) scaled ×0.2 / ×0.7 / ×1.2 |
| 2026-07-14 | `bkgd_pAngle_reweight/` | First version of `pAngle` (fit + toys, no xsec) |
| 2026-09-02 → 09-15 | `test_fit*.root`, `test_xsec*.root` | Tests on `gundam_v27p2.root`, PPFX weight added |

### v30 inputs (`gundam_v30.root`), no sideband, no G4 — `data/<var>/<tag>/`

| Date | Var | Tag | Mean xsec unc | Notes |
|---|---|---|---|---|
| 2026-09-29 | EleE | pAngle | 25.4% | |
| 2026-10-06 | EleE | emb, gainLow, gainHigh | 24.9 / 23.2 / 23.3% | Detector variations from `merged_v30.root` |
| 2026-10-08 → 10-09 | EleE | frac_uncertainty_asimov (+ `_diagonal`, `_noflux`, `_diagonal_noflux`) | — | Systematic-subset breakdown |
| 2026-10-09 | EleE | agent_asimov | 23.1% | Test of the agentic pipeline; matches `asimov` |
| 2026-10-14 | EleE | asimov | 23.1% | |
| 2026-10-14 → 10-16 | EleE | decompose_asimov, split_detector_asimov[_nomcthrow], split_genie_asimov | — | Found no single dominant driver (GENIE ~12, detector ~12, bkg-stat ~10, flux-norm ~9, flux-beam ~8%) |
| 2026-10-16 → 10-19 | EleE | Q2, res1p5, nc, res0p5 | 25.0 / 22.6 / 23.1 / 25.1% | |
| 2026-10-23 | EleE | data_10pct | 34.5% | 10% data |
| 2026-10-23 | EleAngle | asimov | 20.7% | First EleAngle run |

### v31 inputs (`gundam_v31.root`) — `data/<var>/v31/<tag>/`

Batch run on 2026-10-23 11:37–13:12 for all tags and both variables (logs in
`data/v31_logs/`). **No NC π0 sideband, no G4.** `summary.txt` marks the
frac_uncertainty and decompose runs `FAILED(rc=1)`, but their outputs and
plots were all written, so the failure came after the plots.

| Date | Var | Tag | Sideband | G4 | Mean xsec unc | Notes |
|---|---|---|---|---|---|---|
| 10-23 | EleE | nc / Q2 / res0p5 / res1p5 | no | no | 23.2 / 24.7 / 25.0 / 22.6% | |
| 10-23 | EleE | pAngle / emb / gainHigh / gainLow | no | no | 24.9 / 24.2 / 23.0 / 23.9% | Detector variations from `merged_final_v31e_flux.root` |
| 10-23 | EleE | frac_uncertainty ×4, decompose, split_detector ×2, split_genie (+rerun 13:59) | no | no | — | Decompose: quadrature sum 24.1% vs total 23.5% |
| 10-23 | EleAngle | nc / Q2 / res0p5 / res1p5 | no | no | 20.7 / 21.9 / 21.2 / 19.6% | |
| 10-23 | EleAngle | pAngle / emb / gainHigh / gainLow | no | no | 20.6 / 20.9 / 21.1 / 20.3% | |
| 10-23 | EleAngle | asimov | no | no | 20.4% | |
| 10-23 | EleAngle | frac_uncertainty ×4, decompose, split_* | no | no | — | |
| 10-23 16:47 | EleE | asimov | no (`NO_SIDEBAND=1`) | no | 23.5% | Re-run of `asimov`; reference for the sideband comparison |
| 10-27 09:30 | EleAngle | asimov_sb | **yes** | no | 19.6% | NC π0 sideband added (vs 20.4%) |
| 10-27 09:39 | EleE | asimov_sb | **yes** | no | 20.2% | NC π0 sideband added (vs 23.5%) |
| 10-27 10:27 | EleE / EleAngle | data_10pct_sb | yes | no | 32.9 / 33.1% | 10% data with sideband |
| 10-27 18:14 | EleE | asimov_sb_g4 | yes | **yes** | 21.0% | First run with G4 (vs 20.2% for `asimov_sb`). Most of that 0.8-point gap is toy noise (the same config gives 20.4% in the N−1 study); the N−1 study puts G4 at ≈1–2%. Fit status OK, covariance accurate. Per-bin: 25.5, 20.3, 21.5, 20.5, 21.5, 19.8, 20.7, 19.9, 20.6, 19.3%. `OUTPUT=data/EleE/v31/asimov_sb_g4` |
| 10-27 18:20 | EleE | frac_uncertainty_asimov_sb_g4 | yes | yes | 20.4% (all) | N−1 study with the G4 variants (`g4_only`, `no_g4`). Mean contribution alone / within total (N−1): GENIE 11.1 / 11.3%, detector 7.2 / 8.2%, flux 3.8 / 5.9%, **G4 0.8 / 1.7%**; stat-only 13.6%. Removing G4 gives 20.5%, the same as the total within toy noise (1000 toys ≈ ±0.5 points). G4 shows up only in a few bins (N−1 up to 5–7% at 1.25–2 GeV) and is negligible elsewhere. `OUTPUT=data/EleE/v31/frac_uncertainty_asimov_sb_g4` |

### Config changes since the last run (not yet in any study)

- **2026-10-27 — G4 reinteraction systematics added** to
  `parameterSetList_EleE.yaml` (proton, π+, π−). Normalizations binned in true
  leading-pion KE (0/90/210/420 MeV) and leading-proton p (0/400/600/850
  MeV/c), with fractional covariances from 2000 Geant4 universes
  (`Configs_ParameterSet/G4/`). `parameterLimits` removed from the G4 sets (they
  conflicted with eigen-decomposition). Every EleE run from now on includes G4;
  EleAngle doesn't yet.
- **2026-10-27 — G4 added to the N−1 study (EleE).** The
  `frac_uncertainty/parameterSetList_{all,no_xsec,no_flux,no_detector}.yaml`
  lists now include the three G4 sets, and new `g4_only` / `no_g4` lists were
  added. `run_frac_uncertainty.sh` runs those two variants whenever the
  `g4_only` list exists, and `plot_frac_uncertainty.py` draws a G4 curve when
  they're present. EleAngle lists are unchanged. Earlier frac-uncertainty
  results had no G4 in "all".
- **2026-10-27 — `makePrettyPlots.py`** now loops over every parameter set in
  the fit (G4 included) and writes `<set>Summary.png`: pulls on top, covariance
  and correlation below. Existing per-plot file names are unchanged.
- **2026-10-28 — `makePrettyPlots.py` plot changes.** Dropped the
  `<set>Summary.png` figures. The three G4 sets are now drawn together as
  `g4FitPulls.png`, `g4CovarianceMatrix.png` and `g4CorrelationMatrix.png`,
  replacing the per-set G4 plots. The matrices are block-diagonal, from each
  set's own post-fit matrix: the fit stores G4 only eigen-decomposed in the total
  covariance, so there are no physical-basis cross-set terms. G4 bins are labelled
  from each set's binning file (proton momentum 0/400/600/850 MeV/c; π± kinetic
  energy 0/90/210/420 MeV). In the total covariance plot the G4 eigen entries
  read e.g. "G4 Proton eigen 0". Regenerated `plots/EleE/v31/asimov_sb_g4`.

## v31p2 — all fake-data tags with the NC π0 sideband (2026-10-28)

Run 11:43–12:36 with `scratchpad/run_v31p2.sh`. Outputs are in
`data/<var>/v31p2/<tag>` and `plots/<var>/v31p2/<tag>`; logs and `summary.txt`
are in `data/v31p2_logs/`. Both variables use the sideband, with fit/toy
sideband bins of 10/5 (edges 0, 0.33, 0.48, 0.67, 0.96, 5 GeV). EleE includes
G4 and EleAngle doesn't. Every completed fit had status OK and an accurate
covariance. Mean per-bin xsec uncertainty, v31 (no sideband) → v31p2:

| Tag | EleE | EleAngle |
|---|---|---|
| asimov | 23.5 → 20.6% | 20.4 → 19.9% |
| emb | 24.2 → 20.6% | 20.9 → 20.1% |
| gainHigh | 23.0 → 20.4% | 21.1 → 20.3% |
| gainLow | 23.9 → 21.2% | 20.3 → 19.7% |
| nc | 23.2 → 19.9% | 20.7 → 20.0% |
| Q2 | 24.7 → 21.3% | 21.9 → 20.6% |
| res0p5 | 25.0 → 20.5% | 21.2 → 19.9% |
| res1p5 | 22.6 → 20.7% | 19.6 → 18.8% |
| pAngle | 24.9 → 20.6% | 20.6 → 19.1% |

### Config changes for v31p2

- **2026-10-28 — Sideband added to the EleAngle fake-data configs** for emb,
  gainHigh, gainLow, nc, pAngle, Q2, res0p5 and res1p5 (Fitter and Toy). Each
  block is copied verbatim from the matching EleE config, so per-tag fake-data
  reweights (e.g. `nc`) also apply to the angle sideband.
- **2026-10-28 — pAngle sideband weight guarded against NaN.** The pAngle
  fake-data weight p·sin(acos(cosθ_e)) is NaN for the 11503/15194 sideband events
  that have no electron candidate. The sideband data formula now applies the
  weight only where both electron variables are defined, and uses 1 otherwise.
  The user chose this over leaving the sideband unweighted. The sideband events
  that do have an electron candidate have p = 0.3–524 MeV, so their weight is
  about 0.1–0.5.
- **2026-10-28 12:40–12:45 — v31p2 pAngle re-run** with the guarded weight.
  Both fits OK with an accurate covariance. The weight pulls the sideband fake
  data 24% below MC (3875 vs 5096 events), against 9% below for the signal
  (3698 vs 4042). Results: EleE 20.6%, EleAngle 19.1%.


## 2026-10-28 — EleAngle: G4 enabled, v31p2 re-run + N-1 study
- Config: uncommented the three G4 sets (Proton, PiPlus, PiMinus) in `Configs_ParameterSetList/parameterSetList_EleAngle.yaml`, so EleAngle now matches EleE.
- Config: rebuilt `frac_uncertainty/EleAngle/parameterSetList_*.yaml` from the EleE lists (template swapped to TrueEleAngle). Added `g4_only` and `no_g4`; G4 is now in all / no_xsec / no_flux / no_detector.
- The previous EleAngle v31p2 results (sideband, no G4) moved to `data/EleAngle/v31p2_noG4` and `plots/EleAngle/v31p2_noG4`.
- Running: all 9 tags → `data/EleAngle/v31p2/<tag>` (logs `data/v31p2_logs/EleAngle_<tag>_g4.log`), then N-1 → `data/EleAngle/v31/frac_uncertainty_asimov_sb_g4` (plots `plots/EleAngle/v31/frac_uncertainty_asimov_sb_g4`).

## 2026-10-28 — Closure χ² table
- New `scripts/closure_chi2_table.py`: scans `data/<var>/<version>/<study>/xsec_<tag>.root` and tabulates the GUNDAM-vs-CAF-truth correlated and uncorrelated χ²/ndof and p-values. Output: `notes/closure_chi2_table.{md,csv}`. Filter with `--vars` and `--versions`.
- `makePrettyPlots.py`: the χ² calculation moved into `closure_chi2()`, which the comparison plot and the table both call. Numbers are unchanged (EleE v31p2 pAngle: 8.3/10, p = 0.600).
- No closure is bad: the largest is EleAngle emb (v31p2_noG4 19.7/13, p = 0.10). All EleE are χ²/ndof < 0.9.
- `scripts/closure_chi2_table.py` now also writes `plots/closure_chi2_table_{EleE,EleAngle}.png` (`--png-dir` to change), with cells shaded red where p < 0.05.
- Changed: one table per variable and version, `plots/<var>/<version>/closure_chi2_table.{png,md,csv}`. Removed the combined `plots/closure_chi2_table_*.png` and `notes/closure_chi2_table.{md,csv}`. Option is now `--plots-dir`.

## 2026-10-28 — Why the EleE pAngle closure p-value dropped (v31 p=0.93 → v31p2 p=0.60)
- Diagnostic run `data/EleE/v31p2_diag/pAngle_sbNominal`: v31p2 pAngle, but the sideband fake data is nominal MC (pAngle weight only in the signal region). Configs are in the scratchpad (FIT_CONFIG/TOY_CONFIG override). Fit OK, accurate covariance. Closure χ² 10.4/10, p=0.41.
- Closure bias by bin, (GUNDAM − truth)/truth %:
  - v31, no sideband: [−43 −11 −11 −6 −5 −2 −4 −2 −2 +12], integrated −7.9%.
  - Sideband nominal: [−45 −13 −11 −7 −5 −2 −5 −2 −3 +11], integrated −8.6%, the same as v31.
  - v31p2, sideband reweighted (24% low): [−24 +2 −3 +19 +4 +5 +22 +18 +46 +49], integrated +12.8%.
- Cause: mostly the uncertainty, not a new bias. The pAngle bias (lowest bin −45%) exists with or without the sideband. The sideband and G4 shrink the xsec uncertainty, most in the lowest bin (25%→18%).
  - The v31 residual against the sideband-nominal covariance gives χ² 9.7 (4.4 with its own covariance).
- The reweighted sideband adds a second effect. About 20 nuisances get pulled 0.5–1.4σ to absorb the 24% sideband deficit: recombination +1.4σ, NC/CC 2π non-res about −0.8σ, G4 π+ about −1σ. This drags the high-energy templates up (top bin template 1.49) and swaps the low bias for a high one; χ² 8.3.

## 2026-10-28 — EleAngle v31p2 with G4: results
All 9 tags: fit OK, accurate covariance. Xsec uncertainty, mean per bin / integrated %, without G4 (v31p2_noG4) → with G4 (v31p2):
| Tag | no G4 | G4 |
|---|---|---|
| asimov | 19.9/13.7 | 19.5/13.2 |
| emb | 20.1/13.4 | 20.0/13.4 |
| gainHigh | 20.3/14.1 | 19.8/13.4 |
| gainLow | 19.7/13.1 | 19.5/12.8 |
| nc | 20.0/14.0 | 19.3/13.1 |
| pAngle | 19.1/12.7 | 18.6/11.9 |
| Q2 | 20.6/14.0 | 20.6/13.8 |
| res0p5 | 19.9/13.1 | 20.4/13.9 |
| res1p5 | 18.8/12.8 | 19.2/13.2 |
G4 changes these by ±0.5, which is within toy-throw noise. Closure (`plots/EleAngle/v31p2/closure_chi2_table.png`): emb is the worst, 22.0/13 with p = 0.056 (0.103 without G4); everything else has p ≥ 0.70.

## 2026-10-28 — EleAngle N-1 study (sideband + G4)
`data/EleAngle/v31/frac_uncertainty_asimov_sb_g4`, plots `plots/EleAngle/v31/frac_uncertainty_asimov_sb_g4`. All 10 variants OK. Mean/integrated %:
- all 19.9/13.7, template_only 13.8/3.2
- xsec_only 17.8/10.9, flux_only 14.7/5.8, detector_only 14.2/4.6, g4_only 13.8/3.4
- no_xsec 15.6/7.6, no_flux 18.6/11.9, no_detector 19.3/13.1, no_g4 20.1/14.0
GENIE xsec dominates: removing it drops 19.9 → 15.6 (integrated 13.7 → 7.6). Flux and detector are small. G4 is negligible, as for EleE.

## 2026-10-28 — New fake studies: QE / MEC / DIS ×0.5 and ×1.5
- New tags QE0p5, QE1p5, MEC0p5, MEC1p5, DIS0p5, DIS1p5 (true_interaction_mode 0 / 10 / 2), copied from res0p5 with only the weight changed. 24 configs: `fakeStudyConfigs/config_{Fitter,Toy}_FakeData_{EleE,EleAngle}_<tag>.yaml`.
  - The weight applies to the signal and NC π⁰ sideband fake data, not cosmics.
  - The DIS variation also moves the sideband (4,273/15,194 sideband events are DIS; 10,330 are RES).
- `makePrettyPlots.py`: added TRUTH_REWEIGHTS and legend labels for the new tags, so the closure reference uses the reweighted truth.
- `scripts/closure_chi2_table.py` TAG_ORDER, the run-fake-study and run-gundam-fit skill docs, and the tests/run_tests.sh variant lists are updated.
- Running: both variables → `data/<var>/v31p2/<tag>` with sideband + G4. Logs in `data/v31p2_logs/<var>_<tag>.log`.

## 2026-10-28 — QE / MEC / DIS ×0.5 / ×1.5 results (v31p2, sideband + G4)
All 12 runs: fit OK, accurate covariance. Closure χ² vs reweighted truth has p = 1.000 for all of them. Mean per-bin xsec unc %, overall bias of the extracted xsec vs truth:
| Tag | EleE unc | EleE bias | EleAngle unc | EleAngle bias |
|---|---|---|---|---|
| QE0p5 | 21.3 | −1.4% | 20.5 | −2.6% |
| QE1p5 | 20.1 | +0.8% | 18.9 | +1.8% |
| MEC0p5 | 21.4 | 0.0% | 20.6 | +0.3% |
| MEC1p5 | 20.8 | −0.1% | 19.6 | −0.4% |
| DIS0p5 | 20.4 | +8.0% | 18.9 | +1.1% |
| DIS1p5 | 20.8 | −5.9% | 20.3 | +0.5% |
- EleE DIS: a flat normalization bias opposite to the injection, while the EleE signal templates stay close to 1.0. The extracted xsec under-follows the DIS change. Small χ² because it sits inside the correlated (about 13% integrated) uncertainty.
- The sideband (about 28% DIS) pulls the same nuisances by the same amount in both variables: recombination, G4 π+, NC 2π non-res, intrinsic_high, roughly 0.4–0.7σ. Only EleE shows the signal bias; not yet understood.
- Closure tables: `plots/{EleE,EleAngle}/v31p2/closure_chi2_table.png`.

## 2026-10-29 — Closure tables regenerated (v31p2)
- Checked all 30 `data/{EleE,EleAngle}/v31p2/*/fit_*.root`: every one has the G4 sets and the NC π⁰ sideband.
- Re-ran `scripts/closure_chi2_table.py --versions v31p2`; all 15 tags per variable. Numbers unchanged from 10-28. Worst: EleAngle emb 22.0/13 (p = 0.056), EleE pAngle 8.3/10 (p = 0.60); all others p ≥ 0.70.

## 2026-10-29 — G4 universe distributions added
- `Configs_ParameterSet/G4/calc_geant4_covariance.py` now also writes the per-bin universe yield/CV distributions: TH1Ds `<job>/<job>_univ_bin<i>` and `<job>_cv_yield` in the output ROOT file, and `plots/UnivDist_<job>.{pdf,png}`.
- Re-ran all three configs. The covariance matrices are bit-identical to the 10-27 ones, so no fit needs re-running.
- Fractional σ per bin: proton 0.2/2.7/3.0/3.4%; π+ 0.4/1.9/4.5/12.3%; π− 0.1/0.8/3.6/6.3%.

## 2026-10-29 — G4 covariance plots restyled
- New `Configs_ParameterSet/G4/plot_g4_covariance.py`; `calc_geant4_covariance.py` calls it instead of its own plots. The output ROOT now also stores the raw universe yields (`<job>_univ_yields`, TMatrixD) and the binning path.
- Plots in `Configs_ParameterSet/G4/plots/`, png+pdf, in the makePrettyPlots style:
  - `CovMat_<job>` / `CorrMat_<job>`: ROOT kBird / kThermometer with the ICARUS WIP watermark. The axes are titled with the binned variable and the bins are labelled by their edges only.
  - `UnivSpread_<job>`: universe yield/CV vs the binned variable, with all 2000 universes, ±1σ/±2σ bands and the σ per bin.
  - `UnivDist_<job>_bin<i>`: per-bin 1D distribution of universe yield/CV (spineplot nue_icarus style).
- The covariance matrices are unchanged.

## 2026-10-29 — G4 CV vs universe distribution plots
- Replaced `UnivSpread_<job>` with `Configs_ParameterSet/G4/plots/G4Dist_<job>.{png,pdf}`.
  - Top panel: the fine-binned distribution of the binned variable, with the nominal (CV, ppfx_cv_weight) histogram and the ±1σ band of the 2000 G4 universes. Bottom panel: universe/CV. GUNDAM parameter-bin edges are dashed.
  - Bins are 30 MeV for pion KE (0–990) and 50 MeV/c for proton p (0–1500). The last bin includes overflow.
- `plot_g4_covariance.plot_distribution()` is called from `calc_geant4_covariance.py`, or run standalone with `plot_g4_covariance.py <root> --job <job> --config <json>`.
- Check: the CV summed per parameter bin equals the covariance input yields exactly. The covariance matrices are unchanged.
- 380 of 3343 selected events have no leading proton and 2544 of 3343 have no leading pion. Those have no G4 dial and are not plotted.
- The π+ and π− CV histograms are identical: both bin the leading pion KE, whatever its charge. Only the universe band differs.
