function R = evaluateSchedules(schedules, pop, opts)
%EVALUATESCHEDULES Compare imaging or sampling schedules and estimators on a virtual population.
%   R = rlt.evaluateSchedules(schedules, pop)
%
%   schedules : cell array; each cell is a vector of sample times (h)
%   pop       : struct from rlt.simulatePopulation
%
%   For every schedule the same noisy measurements are given to every
%   estimator (a paired comparison). Error is estimate/truth - 1, in percent.
%   For populations without an uptake phase (e.g. the "blood" preset) the hybrid
%   estimator closes the curve before the first sample with Origin="decay" and
%   "select" chooses among mono/biexp/triexp; otherwise Origin="zero" and
%   mono/biexp/uptake.
%
%   Name-value options
%     Methods    estimators from rlt.estimate (default mono, hybrid, stp, select)
%     CV         per-point measurement CV (default 0.10)
%     LambdaPop  population rate for "stp" (default pop.lambdaPop). For an
%                honest result, estimate it on a different cohort (different
%                Seed) than the one being evaluated; see the demo.
%     Seed       RNG seed for the measurement noise
%
%   R is a table with one row per schedule x method:
%     Schedule, Method, NScans, BiasPct (median error), SpreadPct (IQR of error),
%     P90AbsErrPct (90th percentile of |error|), Within10Pct (% of patients
%     with |error| <= 10), NFailed. A failed fit counts as outside 10%.
%     Estimators that need at least two samples ("mono", "select", ...) fail by
%     design on one-sample schedules and show NFailed = pop.N there.
arguments
    schedules cell
    pop (1,1) struct
    opts.Methods (1,:) string = ["mono","hybrid","stp","select"]
    opts.CV (1,1) double {mustBePositive} = 0.10
    opts.LambdaPop (1,1) double = NaN
    opts.Seed (1,1) double = 7
end
if isnan(opts.LambdaPop)
    opts.LambdaPop = pop.lambdaPop;
end
if pop.hasUptake
    origin = "zero";
    models = ["mono","biexp","uptake"];
else
    origin = "decay";
    models = ["mono","biexp","triexp"];
end

sigma  = sqrt(log(1 + opts.CV^2));
stream = RandStream('mt19937ar', 'Seed', opts.Seed);
q      = @rlt.internal.quantileLinear;

nRows    = numel(schedules) * numel(opts.Methods);
Schedule = strings(nRows, 1);
Method   = strings(nRows, 1);
NScans   = zeros(nRows, 1);
BiasPct = zeros(nRows, 1);
SpreadPct = zeros(nRows, 1);
P90AbsErrPct = zeros(nRows, 1);
Within10Pct = zeros(nRows, 1);
NFailed = zeros(nRows, 1);

row = 0;
for s = 1:numel(schedules)
    ts    = schedules{s}(:);
    truth = pop.eval(ts);                                  % N x m
    noisy = truth .* exp(sigma * randn(stream, size(truth)) - sigma^2/2);
    label = strjoin(string(ts'), ",") + " h";

    for m = 1:numel(opts.Methods)
        est = nan(pop.N, 1);
        for i = 1:pop.N
            try
                est(i) = rlt.estimate(opts.Methods(m), ts, noisy(i, :)', ...
                    Isotope=pop.isotope, LambdaPop=opts.LambdaPop, ...
                    Origin=origin, Models=models);
            catch
                est(i) = NaN;
            end
        end
        err = 100 * (est ./ pop.trueTIA - 1);

        row = row + 1;
        Schedule(row)     = label;
        Method(row)       = opts.Methods(m);
        NScans(row)       = numel(ts);
        BiasPct(row)      = q(err, 0.5);
        SpreadPct(row)    = q(err, 0.75) - q(err, 0.25);
        P90AbsErrPct(row) = q(abs(err), 0.9);
        Within10Pct(row)  = 100 * mean(abs(err) <= 10);
        NFailed(row)      = sum(isnan(err));
    end
end

R = table(Schedule, Method, NScans, BiasPct, SpreadPct, P90AbsErrPct, Within10Pct, NFailed);
end
