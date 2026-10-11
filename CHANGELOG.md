# Changelog

## Unreleased
- Released under the MIT licence.
- First run on real MATLAB (CI, R2026b) found and fixed three import bugs: a validation table column named `Row` (reserved by MATLAB tables), CSV column types guessed wrongly by `readtable` (now declared in `rlt.internal.readCsvTable`), and the delimiter guessed wrongly for a file with long quoted text (now fixed to comma). Unit tests: 56 pass.

## 0.2.0
- **Bone-marrow module** (`rlt.marrow.*`): RMBLR, blood-based and image-based self-dose with per-model comparison and optional time-dependent ratio, remainder of body, S-value scaling, total-dose assembly, Monte Carlo uncertainty with shared calibration error, elasticity.
- **Data layer** (`rlt.io.*`): validated import of activity tables (units, mandatory decay-correction flag, duplicates, unit mixing), patient covariates, and S-values with mandatory provenance. Templates in `data/templates/`.
- `rlt.fitTAC`: new `triexp` model; quality diagnostics (extrapolated fraction of the TIA, bound hits, exact determination, poor fit) returned in `fit.warnings`.
- `rlt.tiaHybrid`: `Origin="decay"` for blood and whole-body curves that start at their maximum.
- `rlt.simulatePopulation`: general K-component cohorts; new `blood` preset. **Behaviour change:** the kidney and tumour presets now draw component fractions multiplicatively, so a given seed no longer reproduces v0.1 populations; fields `a1, a2, lambda1, lambda2` became `a`, `lambda` (N x K).
- `rlt.evaluateSchedules` handles populations without uptake (blood): `Origin="decay"` and mono/biexp/triexp selection.
- `rlt.selfCheck`, `rlt.version`, `rlt.provenance`, `rlt.report.marrow`, decay-correction helpers.
- Tests: three suites (`TestRLT`, `TestMarrow`, `TestIO`).
- Docs: `docs/MARROW.md` (methods, assumptions register, literature pitfalls), `docs/DATA_CHECKLIST.md`.

## 0.1.0
- Initial sparse-sampling TIA toolbox: constrained mono/bi-exponential/uptake fits, AICc selection, hybrid and single-time-point estimators, propagation of measurement error, local-deposition dose, virtual cohorts and schedule benchmark.
