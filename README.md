# rlt-sparse-dosimetry

MATLAB toolbox for **time-integrated activity (TIA) and absorbed dose from sparse post-therapy measurements** in radioligand therapy (¹⁷⁷Lu-PSMA, ¹⁷⁷Lu-DOTATATE, ⁹⁰Y, ¹³¹I), with a focus on **red-bone-marrow dosimetry**.

## The problem

Full dosimetry needs many blood samples and scans per cycle. Most centres can afford a few. With so few points, the reported dose depends on choices that rarely get reported: which kinetic model is fitted, how the curve is closed before the first and after the last sample, which time points were acquired, how noisy each point is, and, for marrow, whether blood is a valid proxy at all.

Published work illustrates the stakes: the same ¹⁷⁷Lu-DOTATATE patients gave very different blood-based marrow doses with mono-, bi- and tri-exponential fits, and blood-based doses were well below image-based ones ([details and references](docs/MARROW.md)). This toolbox makes those choices explicit, measurable and auditable.

## What it does

1. **Fits physically constrained models** (mono, bi, tri-exponential, uptake) so effective half-life never exceeds the physical one and activity never goes negative, and flags what the fit cannot tell you (how much of the TIA is extrapolated, exactly determined models, no resolvable washout).
2. **Refuses to over-fit.** AICc selection lets a model compete only when the data can support it.
3. **Marrow dosimetry** (`rlt.marrow`): blood-based and image-based self-dose, a dose for every candidate blood model, time-dependent marrow-to-blood ratio, remainder of body, assembly of all contributions with S-values you supply, and Monte Carlo uncertainty with shared calibration error.
4. **Benchmarks sampling schedules** on virtual cohorts with known truth: for each schedule and estimator, how often is the TIA within 10%?
5. **Validated data intake** (`rlt.io`): units, a mandatory decay-correction flag, duplicates, haematocrit entered as percent, S-values without a recorded reference, all rejected with a complete list of problems.
6. **Audit trail:** every result carries its assumptions, warnings and software provenance; `rlt.report.marrow` writes them to a Markdown report.

## Quick start

Requires MATLAB R2021a or newer. Base MATLAB only, no toolboxes.

```matlab
addpath(pwd)                                     % from the repo root
rlt.selfCheck                                    % known-answer checks, no framework needed
results = runtests("tests")                      % full unit tests
run("examples/demo_marrow_workflow.m")           % end-to-end marrow dosimetry on synthetic data
run("examples/demo_marrow_blood_schedules.m")    % which blood sampling plan is enough?
run("examples/demo_single_patient.m")            % organ/lesion TIA with uncertainty
run("examples/demo_schedule_comparison.m")       % scan schedules for organs and lesions
```

```matlab
% Blood-based marrow self-dose
T   = rlt.io.readActivityTable("activity.csv");              % validated, decay correction undone
P   = rlt.io.readPatientInfo("patients.csv");
b   = rlt.io.getSeries(T, "P001", 1, "blood");
res = rlt.marrow.selfDose(b.t, b.A, ConcentrationUnit=b.unit, Hct=0.45, Sex="M", ...
                          AdministeredMBq=7400, Isotope=b.isotope);
res.doseGy, res.candidates, res.warnings, res.assumptions
rlt.report.marrow("P001_marrow.md", res, Patient="P001");
```

## Functions

| Area | Function | Purpose |
|---|---|---|
| Fitting | `rlt.fitTAC` | constrained model fit, TIA, AICc, quality warnings |
| | `rlt.selectModel` | AICc selection with sparse-data fallback |
| | `rlt.tiaHybrid`, `rlt.tiaSinglePoint`, `rlt.optimalSinglePointTime` | non-fit estimators |
| | `rlt.estimate` | one entry point for all estimators |
| Marrow | `rlt.marrow.rmblr`, `referenceMass` | red-marrow-to-blood ratio, reference mass |
| | `rlt.marrow.selfDose` | blood or VOI based self-dose, all-model comparison |
| | `rlt.marrow.remainderOfBody`, `scaleSValue` | remainder-of-body TIA, S-value mass scaling |
| | `rlt.marrow.totalDose`, `uncertainty`, `elasticity` | sum of contributions, Monte Carlo, sensitivity |
| Data | `rlt.io.readActivityTable`, `validateActivityTable`, `getSeries` | import and validation |
| | `rlt.io.readPatientInfo`, `readSValues`, `lookupS` | covariates and S-values |
| | `rlt.undoDecayCorrection`, `applyDecayCorrection` | decay-correction helpers |
| Dose | `rlt.absorbedDose` | local-deposition dose, or your S-value |
| Study design | `rlt.simulatePopulation`, `evaluateSchedules`, `plotSchedules` | virtual cohorts (kidney, tumour, blood) and schedule benchmark |
| Uncertainty | `rlt.propagateTIA` | measurement error into TIA |
| Housekeeping | `rlt.selfCheck`, `version`, `provenance`, `report.marrow`, `isotope` | |

## Documentation

- [docs/MARROW.md](docs/MARROW.md): methods, assumptions register, when blood-based dosimetry fails, literature.
- [docs/DATA_CHECKLIST.md](docs/DATA_CHECKLIST.md): what to collect so future patient data imports cleanly.
- [docs/DESIGN.md](docs/DESIGN.md): mathematics, design decisions, validation plan, roadmap.
- [CHANGELOG.md](CHANGELOG.md).

## Read this before using it on patients

- **Research and teaching code. Not a medical device and not validated for clinical dosimetry.** It is released under the [MIT licence](LICENSE), which provides it "as is", without warranty. Validating it for any clinical or regulatory use is your responsibility.
- The virtual-cohort presets (`kidney`, `tumour`, `blood`) are **illustrative numbers, not measured data**. Schedule rankings are only as good as the cohort you simulate.
- The toolbox ships **no S-values**. Cross-dose terms need values from OLINDA/EXM, MIRDcalc, IDAC-Dose or your own Monte Carlo, with the reference recorded.
- The blood-based method assumes no specific marrow uptake and no marrow involvement. See [docs/MARROW.md](docs/MARROW.md) before applying it to patients with bone metastases or to ¹⁷⁷Lu-DOTATATE.
- Isotope constants and the blood density are reference values that were not independently verified. Check them (assumptions register in MARROW.md).
- The EANM bone-marrow guideline text was not available when the formulas were implemented; compare them with it.
- The code was written without access to MATLAB. Its numerics were cross-checked in Python, but the MATLAB itself has not been run. Start with `rlt.selfCheck` and `runtests("tests")`; the included GitHub Actions workflow runs the tests on real MATLAB on every push.

## Licence and feedback

[MIT](LICENSE), copyright (c) 2026 Sobhan Shahbazi. Reviews, suggestions and corrections to the formulas, constants and assumptions are welcome as [GitHub issues](https://github.com/callmesobhan/rlt-sparse-dosimetry/issues). The assumptions register in [docs/MARROW.md](docs/MARROW.md) lists what most needs a second pair of eyes.
