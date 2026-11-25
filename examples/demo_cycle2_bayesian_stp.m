%% Cycle-1 -> Cycle-2+ Empirical Bayes NLME Sparse-Schedule Optimizer & Blood-Proxy Audit
% Demonstrates how a clinic ingests its own Cycle-1 multi-time-point patients,
% fits an empirical log-normal population-PK prior, optimizes Cycle-2+ single-time-point
% windows with marrow toxicity-miss bounds, and flags patients where blood-based
% marrow dosimetry violates constant-RMBLR assumptions.

addpath(fileparts(fileparts(mfilename('fullpath'))));

iso = rlt.isotope("Lu177");
stream = RandStream('mt19937ar', 'Seed', 2026);
tFull = [4; 24; 72; 168];

% 1. Build a clinic's Cycle-1 multi-time-point cohort (12 patients)
series = cell(12, 1);
for i = 1:12
    lamBio = 0.015 * exp(0.30 * randn(stream));
    a0     = 110 * exp(0.25 * randn(stream));
    series{i} = struct('t', tFull, 'A', a0 * exp(-(iso.lambdaPhys + lamBio) * tFull));
end

prior = rlt.nlme.fitPopPrior(series, Isotope="Lu177");
fprintf("Fitted empirical Cycle-1 prior from %d patients (lambda_pop = %.4f /h)\n", ...
    prior.nPatients, prior.lambdaPop);

% 2. Optimize Cycle-2+ sparse sampling windows with worst-case toxicity-miss risk
opt = rlt.nlme.optimizeCycle2Schedule(prior, {[24], [48], [72], [96], [120], [24 96]}, ...
    NVirtual=120, DosePerTIA=4.5e-4, DoseLimitGy=2.0, Seed=7);
disp(opt.table);
fprintf("Recommended Cycle-2+ schedule: %s (P90 |err| = %.1f%%, tox-miss = %.1f%%)\n\n", ...
    opt.bestSchedule, opt.bestP90AbsErr, opt.bestToxMissPct);

% 3. Audit blood-based vs vertebral-VOI marrow dosimetry for specific marrow uptake
tB = [2; 24; 48; 96; 168];
cB = 2e-3 * exp(-(iso.lambdaPhys + 0.03) * tB);
cM = 2.0 * rlt.marrow.rmblr(0.45) * 2e-3 * exp(-(iso.lambdaPhys + 0.006) * tB);
audit = rlt.marrow.auditBloodProxy(tB, cB, tB, cM, Hct=0.45);
fprintf("Blood-proxy audit: %s (image/blood dose ratio = %.2fx)\n  %s\n", ...
    audit.status, audit.doseRatio, audit.recommendation);
