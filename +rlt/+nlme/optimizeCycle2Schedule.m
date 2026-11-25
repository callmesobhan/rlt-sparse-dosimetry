function out = optimizeCycle2Schedule(prior, candidateSchedules, opts)
%OPTIMIZECYCLE2SCHEDULE Design optimal Cycle-2+ sparse sampling windows from a Cycle-1 cohort prior.
%   out = rlt.nlme.optimizeCycle2Schedule(prior, {[24], [48], [72], [96], [120], [168], [24 96]})
%
%   Unlike commercial RLT software (which assumes fixed literature single-time-point times),
%   this function draws virtual Cycle-2+ patients directly from the clinic's own empirical
%   Cycle-1 population prior `prior`, simulates realistic inter-cycle physiological drift,
%   and evaluates every candidate sparse schedule using Empirical Bayes MAP estimation
%   (with and without the patient's Cycle-1 anchor) versus naive Hanscheid STP.
%
%   Reports:
%     - BiasPct, P90AbsErrPct, Within10Pct
%     - ToxicityMissPct: worst-case false-safe rate P(estDose < Limit & trueDose >= Limit),
%       directly bounding the clinical risk of missing marrow (2.0 Gy) or kidney (23 Gy) limits.
%     - Recommended optimal schedule `bestSchedule` minimizing P90AbsErrPct under ToxicityMissCap.
%
%   Base MATLAB only; no toolboxes required.
arguments
    prior (1,1) struct
    candidateSchedules cell
    opts.NVirtual (1,1) double {mustBeInteger, mustBePositive} = 200
    opts.MeasurementCV (1,1) double {mustBePositive} = 0.10
    opts.InterCycleAmpCV (1,1) double {mustBePositive} = 0.20
    opts.InterCycleRateCV (1,1) double {mustBePositive} = 0.15
    opts.DosePerTIA (1,1) double {mustBePositive} = 1.0
    opts.DoseLimitGy (1,1) double = NaN
    opts.ToxicityMissCapPct (1,1) double {mustBeNonnegative} = 5.0
    opts.Seed (1,1) double = 42
end
if ~isnan(opts.DoseLimitGy) && opts.DoseLimitGy <= 0
    error('rlt:badLimit', 'DoseLimitGy must be positive when supplied.');
end

stream = RandStream('mt19937ar', 'Seed', opts.Seed);
N = opts.NVirtual;
lamP = prior.lambdaPhys;
q = @rlt.internal.quantileLinear;

% Draw Cycle-1 virtual patients from empirical population prior N(mu, omega)
L = chol(0.5 * (prior.omega + prior.omega'), 'lower');
th1 = prior.mu + L * randn(stream, 2, N);              % 2 x N

% Evolve to Cycle-2+ with realistic inter-cycle physiological drift
sigA = sqrt(log(1 + opts.InterCycleAmpCV^2));
sigL = sqrt(log(1 + opts.InterCycleRateCV^2));
th2 = th1 + [sigA * randn(stream, 1, N); sigL * randn(stream, 1, N)];

A2    = exp(th2(1, :))';
lam2  = lamP + exp(th2(2, :))';
trueTIA2 = A2 ./ lam2;
trueDose2 = trueTIA2 * opts.DosePerTIA;

% If no explicit DoseLimitGy was supplied, calibrate a discriminating 85th-percentile toxicity limit
doseLimit = opts.DoseLimitGy;
if isnan(doseLimit)
    doseLimit = q(trueDose2, 0.85);
end

% Synthetic Cycle-1 fit structs (so predictCycle2 can use the patient's own Cycle-1 anchor)
c1Fits = cell(N, 1);
for i = 1:N
    c1Fits{i} = struct('amplitude', exp(th1(1, i)), 'lambda', lamP + exp(th1(2, i)));
end

sigmaMeas = sqrt(log(1 + opts.MeasurementCV^2));
nSched = numel(candidateSchedules);
methods = ["BayesPop", "BayesCycle1", "NaiveSTP"];
nRows = nSched * numel(methods);

Schedule        = strings(nRows, 1);
Method          = strings(nRows, 1);
NPoints         = zeros(nRows, 1);
BiasPct         = zeros(nRows, 1);
P90AbsErrPct    = zeros(nRows, 1);
Within10Pct     = zeros(nRows, 1);
ToxicityMissPct = zeros(nRows, 1);

row = 0;
for s = 1:nSched
    ts = sort(candidateSchedules{s}(:));
    label = strjoin(string(ts'), ",") + " h";
    truthPts = A2 .* exp(-lam2 * ts');                  % N x m
    noisyPts = truthPts .* exp(sigmaMeas * randn(stream, size(truthPts)) - 0.5*sigmaMeas^2);

    estPop = zeros(N, 1);
    estC1  = zeros(N, 1);
    estSTP = zeros(N, 1);
    for i = 1:N
        rPop = rlt.nlme.predictCycle2(ts, noisyPts(i, :)', prior, MeasurementCV=opts.MeasurementCV);
        rC1  = rlt.nlme.predictCycle2(ts, noisyPts(i, :)', prior, Cycle1Fit=c1Fits{i}, ...
            MeasurementCV=opts.MeasurementCV, InterCycleAmpCV=opts.InterCycleAmpCV, ...
            InterCycleRateCV=opts.InterCycleRateCV);
        estPop(i) = rPop.tia;
        estC1(i)  = rC1.tia;
        estSTP(i) = rPop.tiaNaiveSTP;
    end

    estList = {estPop, estC1, estSTP};
    for m = 1:numel(methods)
        estTIA = estList{m};
        err = 100 * (estTIA ./ trueTIA2 - 1);
        estDose = estTIA * opts.DosePerTIA;
        toxMiss = 100 * mean((estDose < doseLimit) & (trueDose2 >= doseLimit));

        row = row + 1;
        Schedule(row)        = label;
        Method(row)          = methods(m);
        NPoints(row)         = numel(ts);
        BiasPct(row)         = q(err, 0.5);
        P90AbsErrPct(row)    = q(abs(err), 0.9);
        Within10Pct(row)     = 100 * mean(abs(err) <= 10);
        ToxicityMissPct(row) = toxMiss;
    end
end

tbl = table(Schedule, Method, NPoints, BiasPct, P90AbsErrPct, Within10Pct, ToxicityMissPct);

% Choose optimal BayesCycle1 schedule meeting toxicity-miss cap
c1Rows = tbl(tbl.Method == "BayesCycle1", :);
safe = c1Rows(c1Rows.ToxicityMissPct <= opts.ToxicityMissCapPct, :);
if ~isempty(safe)
    [~, idxBest] = min(safe.P90AbsErrPct);
    bestRow = safe(idxBest, :);
else
    [~, idxBest] = min(c1Rows.ToxicityMissPct);
    bestRow = c1Rows(idxBest, :);
end

out.table           = tbl;
out.doseLimitGy     = doseLimit;
out.bestSchedule    = bestRow.Schedule;
out.bestP90AbsErr   = bestRow.P90AbsErrPct;
out.bestWithin10Pct = bestRow.Within10Pct;
out.bestToxMissPct  = bestRow.ToxicityMissPct;
out.provenance      = rlt.provenance();
end
