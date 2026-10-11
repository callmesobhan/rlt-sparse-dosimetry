# Design notes (v0.2)

## 1. Problem statement

**Who:** a nuclear-medicine physics group running ¹⁷⁷Lu (or ⁹⁰Y, ¹³¹I) radioligand therapy that can afford only 1–3 quantitative scans per cycle.

**Question:** given scans at times `t₁…tₙ`, what is the time-integrated activity (TIA) of an organ or lesion, how uncertain is it, and which scan schedule makes it trustworthy enough to act on?

**Why it is hard:** TIA = ∫₀^∞ A(t) dt. Scans cover only part of the curve, so the integral is dominated by an assumption: the shape before the first scan and after the last. Expect the choice of fit model and time points to move TIA substantially, yet clinical reports quote a single number. The benchmark in this toolbox exists to measure by how much, for your own cohort, instead of assuming.

**Deliverable:** a base-MATLAB toolbox that (a) estimates TIA and dose from sparse data with explicit, physically sensible assumptions, (b) quantifies uncertainty, (c) benchmarks schedules and estimators against known truth on a virtual cohort.

**Primary use case (v0.2):** red-bone-marrow dosimetry for ¹⁷⁷Lu radioligand therapy, where the quantity of interest cannot be imaged directly. See [MARROW.md](MARROW.md).

**Out of scope:** image processing, registration, segmentation, voxel dosimetry, Monte Carlo transport, spongiosa micro-structure models, α-emitters with decay chains (²²⁵Ac).

## 2. Mathematics

Let λp be the physical decay constant. Every fitted rate is `λp + biological term ≥ 0`, so no fit implies an effective half-life longer than the physical one.

| Model | A(t) | TIA | Free parameters |
|---|---|---|---|
| `mono` | a·e^(−λ₁t) | a/λ₁ | 2 |
| `biexp` | a₁e^(−λ₁t) + a₂e^(−λ₂t), λ₁<λ₂ | a₁/λ₁ + a₂/λ₂ | 4 |
| `triexp` | three exponentials, λ₁<λ₂<λ₃ | Σ aₖ/λₖ | 6 |
| `uptake` | a(e^(−λ₁t) − e^(−λ₂t)), λ₁<λ₂ | a(1/λ₁ − 1/λ₂) | 3 |

**Fitting by variable projection.** Amplitudes are linear and must be ≥ 0, so for each trial set of rates they are solved exactly by `lsqnonneg`. Only the rates are searched (`fminbnd` for one rate, `fminsearch` with 6 starts for two, 8 for three). Rates are parametrised as `λ₁ = λp + eˣ`, `λ₂ = λ₁ + eʸ`, which enforces the ordering and bound without constrained optimisation, so no Optimization Toolbox is needed.

**Weights.** Default residual scaling is 1/A (constant relative error), a reasonable first model for SPECT quantification error.

**Model selection.** AICc = n·ln(SSE/n) + 2p + 2p(p+1)/(n−p−1). It is undefined unless n ≥ p+2. So with 3 scans nothing can be compared (mono needs n ≥ 4, uptake n ≥ 5, biexp n ≥ 6). `selectModel` therefore returns mono-exponential whenever no model has a finite AICc. **This is a result, not a limitation:** with typical clinical schedules the data cannot choose a model, which is exactly why the population-prior ideas in the roadmap matter.

**Fit diagnostics.** Each fit reports the share of its TIA that lies outside the sampled window (computed exactly from the fitted exponentials), whether the model is exactly determined, whether the slowest rate sits on the physical-decay bound (no washout resolved), and whether the misfit is large. They are returned in `fit.warnings` and copied into marrow results, so a number never travels without its caveats.

**Single time point.** If true decay is mono-exponential with rate λ and we assume λ̂, the TIA estimated from one scan at time `ts` is `A(ts)·e^(λ̂·ts)/λ̂`. Its relative error is `(λ/λ̂)·e^((λ̂−λ)·ts) − 1`. The derivative with respect to λ̂ at λ̂ = λ is `ts − 1/λ`, which vanishes at `ts = 1/λ = T_eff/ln 2`. For a kidney T_eff of about 50 h this is about 3 days, consistent with the single-time-point literature for ¹⁷⁷Lu.

**Dose (local deposition).** `D[Gy] = TIA[MBq·h] · 3.6×10⁹ · E_β[MeV] · 1.602×10⁻¹³ · φ / m[kg]`. Photons and cross-dose ignored; φ is user-supplied.

**Uncertainty.** `propagateTIA` multiplies each measurement by mean-one lognormal noise (coefficient of variation CV), re-estimates, and reports the 2.5 / 50 / 97.5 percentiles. It captures measurement error for the given schedule only. It starts from already-noisy data, so it is a first-order propagation, not a posterior. Model-form error is assessed separately by `evaluateSchedules`.

## 3. Virtual cohort

`simulatePopulation` generates `A(t) = a₁e^(−λ₁t) + a₂e^(−λ₂t) − (a₁+a₂)e^(−k_u t)`. With `λ₁ ≤ λ₂ ≤ k_u/1.2` this is non-negative and starts at zero. The truth is closed-form: `a₁/λ₁ + a₂/λ₂ − (a₁+a₂)/k_u`.

The generating model (two washout phases plus uptake) is deliberately richer than the estimators' mono-exponential, so that model-form bias shows up in the benchmark. Parameter distributions live in `+rlt/+internal/presetParams.m` and are **illustrative**.

## 4. Layout

```
+rlt/                 public API (package: call as rlt.fitTAC, ...)
+rlt/+marrow/         red-marrow dosimetry (rmblr, selfDose, remainderOfBody, totalDose, ...)
+rlt/+io/             validated data import (activity tables, patient covariates, S-values)
+rlt/+report/         Markdown report with assumptions and provenance
+rlt/+internal/       helpers not meant for direct use
data/templates/       CSV templates with a synthetic example patient
examples/             runnable demos
tests/                TestRLT, TestMarrow, TestIO (matlab.unittest)
docs/                 MARROW.md, DATA_CHECKLIST.md, DESIGN.md
.github/workflows/    CI that runs the tests on real MATLAB
```

## 5. Decisions

| Decision | Reason |
|---|---|
| Package namespace `+rlt` | No name clashes with user code; `rlt.` prefix reads clearly in scripts |
| Base MATLAB only | Clinical physics groups often lack the Optimization and Statistics Toolboxes |
| `arguments` blocks, name=value | Early, readable errors; needs R2021a+ |
| One `rlt.estimate` dispatcher | Benchmarks and uncertainty code treat all estimators identically |
| Private RNG streams | Seeding never disturbs the user's global `rng` state |
| Struct outputs, no classes yet | Simplest thing that works; revisit if the API grows |
| No silent defaults for patient-specific or uncertain inputs (Hct, uncertainty CVs, S-values, decay-correction flag) | A default that is right for one centre is a hidden error for another; missing inputs raise errors that say what to supply |
| Collect all validation problems, then fail once | One round trip to fix a bad data file instead of one error at a time |
| Every result carries `assumptions`, `warnings`, `provenance` | A dose without its assumptions cannot be audited or compared |
| Marrow self-dose via concentration, not mass | Under local deposition the marrow mass cancels, so an unmeasurable patient quantity is not needed |

## 6. Validation plan

1. **Analytic tests** (in `tests/`): exact recovery of mono-, bi- and tri-exponential TIA, closed-form extrapolated fraction, trapezoid convergence, single-point optimum at 1/λ, dose constant, marrow-mass independence, analytic elasticities (Hct/(1−Hct), 1, 1), cohort truth against numerical integration, error paths for every validator.
2. **Cross-implementation check (done for v0.2):** the core numerics were mirrored in Python (NumPy/SciPy) and agreed with the test expectations. This does not replace running MATLAB: use `rlt.selfCheck` and `runtests("tests")`, and the CI workflow.
3. **Property tests (to add):** fitted rates ≥ λp for random data; TIA monotone in tail rate; AICc ordering with noise-free generated data; round trips through the CSV layer for random tables.
4. **Real-data comparison (needs patient data):** for patients with a full schedule, treat the full-schedule result as the reference, subsample, and compare against what `evaluateSchedules` predicts. For marrow, compare blood-based against image-based and against toxicity (platelet nadir), as in the cited studies.
5. **Cross-check against validated software:** compare absorbed dose with OLINDA/EXM, MIRDcalc or your department's dosimetry software on shared cases.
6. **Formula review:** check `selfDose`, `remainderOfBody` and the S-value handling line by line against the EANM bone-marrow guideline (not available when v0.2 was written).

## 7. Roadmap

- **Next, once any patient data exist:** run `docs/DATA_CHECKLIST.md` on the first patient end to end; replace the illustrative presets with distributions fitted to your own cohort; add the dose–response step (marrow dose versus platelet nadir).
- **Population-prior fitting (hierarchical Bayesian):** borrow curve shape from the cohort, individual scale from the scans; this is what makes 2–3 samples informative and gives honest intervals for blood-based dosimetry.
- **Active-marrow and lesion masks as inputs:** per-VOI lesion fractions that correct the image-based marrow activity instead of only warning.
- **S-value helpers:** import from MIRDcalc / IDAC-Dose export formats; mass- and sex-specific scaling rules with documented exponents.
- **Voxel-level TIA maps** with a single time point; lesion sphere S-values; ²²⁵Ac chain.
- **Packaging:** `.mltbx` toolbox file, function signature files for autocomplete, documentation as MATLAB live scripts.

## 8. Open questions for the owner

1. Which agent first: ¹⁷⁷Lu-PSMA-617 or ¹⁷⁷Lu-DOTATATE? They differ in the likelihood of specific marrow uptake and in bone-metastasis prevalence, which decides the default method.
2. Does your imaging protocol include a late whole-body acquisition (needed for remainder of body)?
3. Which S-value source will you use (OLINDA/EXM, MIRDcalc, IDAC-Dose, Monte Carlo)?
4. Licence (none is set; add one before sharing the repository).
5. Is MATLAB R2021a the right minimum for your collaborators?
