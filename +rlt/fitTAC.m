function fit = fitTAC(t, A, opts)
%FITTAC Fit a physically constrained kinetic model to a sparse time-activity curve.
%   fit = rlt.fitTAC(t, A) fits a mono-exponential model.
%   fit = rlt.fitTAC(t, A, Model="biexp", Isotope="Lu177")
%
%   t  : sample times after administration (h), strictly increasing
%   A  : measured activity or concentration, NOT decay-corrected
%        (use rlt.undoDecayCorrection first if your data are decay-corrected)
%
%   Models (lambda_phys = physical decay constant of the isotope):
%     "mono"   A(t) = a*exp(-l1*t)                              l1 >= lambda_phys
%     "biexp"  A(t) = a1*exp(-l1*t) + a2*exp(-l2*t)             lambda_phys <= l1 < l2
%     "triexp" A(t) = sum of three exponentials                 lambda_phys <= l1 < l2 < l3
%     "uptake" A(t) = a*(exp(-l1*t) - exp(-l2*t))               lambda_phys <= l1 < l2
%   Every rate is lambda_phys plus a non-negative biological term, and every
%   amplitude is >= 0, so a fit can never imply an effective half-life longer
%   than the physical one or negative activity. Use "mono"/"biexp"/"triexp"
%   for curves that start high and fall (blood, whole body) and "uptake" for
%   organs and lesions that fill first.
%
%   Name-value options
%     Isotope          "Lu177" (default), "Y90", "I131"
%     Model            "mono" (default), "biexp", "triexp", "uptake"
%     Weights          "relative" (default, residuals scaled by 1/A) or "none"
%     WarnExtrapolated fraction of TIA lying outside [t(1), t(end)] above which a
%                      warning is recorded (default 0.5, a heuristic)
%
%   Output struct fields
%     model, isotope, n, nParams, weights, t, A
%     amplitude   fitted amplitudes (>= 0)
%     lambda      fitted rates (1/h)
%     tia         time-integrated activity, 0..Inf (units of A times h)
%     sse         weighted sum of squared residuals
%     weightedRMSE  sqrt(sse/n); with relative weights, the typical fractional misfit
%     aicc        corrected AIC; NaN when n < nParams + 2
%     exactlyDetermined  true when n == nParams (zero residual, no df)
%     extrapolatedFraction  share of TIA outside the measured window [t(1), t(end)]
%     warnings    string array of quality warnings (empty when none)
%     predict     function handle: predict(tq) -> fitted activity at times tq
%
%   Amplitudes are solved by non-negative least squares for each trial set of
%   rates (variable projection), so only the rates are searched numerically.
%   Base MATLAB only: no toolboxes required.
arguments
    t double {mustBeVector, mustBeNonnegative, mustBeFinite}
    A double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Isotope (1,1) string = "Lu177"
    opts.Model (1,1) string {mustBeMember(opts.Model, ["mono","biexp","triexp","uptake"])} = "mono"
    opts.Weights (1,1) string {mustBeMember(opts.Weights, ["relative","none"])} = "relative"
    opts.WarnExtrapolated (1,1) double = 0.5
end
t = t(:);
A = A(:);
rlt.internal.checkSeries(t, A);

iso   = rlt.isotope(opts.Isotope);
lamP  = iso.lambdaPhys;
model = opts.Model;
n     = numel(t);
p     = rlt.internal.nParams(model);
if n < p
    error('rlt:tooFewPoints', ...
        'Model "%s" has %d parameters but only %d time points were given.', model, p, n);
end

if opts.Weights == "relative"
    wts = 1 ./ max(A, 0.01*max(A));
else
    wts = ones(n, 1);
end

obj = @(th) projectedSSE(model, th, t, A, wts, lamP);
th  = searchRates(model, obj);
[sse, a, rates] = obj(th);

fit.model     = model;
fit.isotope   = iso.name;
fit.n         = n;
fit.nParams   = p;
fit.weights   = opts.Weights;
fit.t         = t;
fit.A         = A;
fit.amplitude = a;
fit.lambda    = rates(:);
fit.tia       = tiaFromParams(model, a, rates);
fit.sse       = sse;
fit.weightedRMSE = sqrt(sse / n);
if n - p - 1 > 0
    floorSSE = 1e-12 * sum((wts .* A).^2);
    fit.aicc = n*log(max(sse, floorSSE)/n) + 2*p + 2*p*(p+1)/(n-p-1);
else
    fit.aicc = NaN;
end
fit.exactlyDetermined = (n == p);

window = windowIntegral(model, a, rates, t(1), t(end));
if fit.tia > 0
    fit.extrapolatedFraction = min(max(1 - window/fit.tia, 0), 1);
else
    fit.extrapolatedFraction = NaN;
end
fit.warnings = qualityWarnings(fit, lamP, opts.WarnExtrapolated);
fit.predict  = @(tq) basisMatrix(model, rates, tq(:)) * a;
end

% ---------------------------------------------------------------------------

function [sse, a, rates] = projectedSSE(model, th, t, A, w, lamP)
rates = thetaToRates(model, th, lamP);
X = basisMatrix(model, rates, t);
[a, sse] = lsqnonneg(w .* X, w .* A);
end

function rates = thetaToRates(model, th, lamP)
th = min(max(th, log(1e-6)), log(20));     % keep rates in a sane range
l1 = lamP + exp(th(1));                    % biological term >= 0
switch model
    case "mono"
        rates = l1;
    case {"biexp","uptake"}
        rates = [l1, l1 + exp(th(2))];     % l2 > l1 by construction
    case "triexp"
        l2 = l1 + exp(th(2));
        rates = [l1, l2, l2 + exp(th(3))];
end
end

function X = basisMatrix(model, rates, t)
switch model
    case "mono"
        X = exp(-t * rates(1));
    case {"biexp","triexp"}
        X = exp(-t * rates);
    case "uptake"
        X = exp(-t * rates(1)) - exp(-t * rates(2));
end
end

function tia = tiaFromParams(model, a, rates)
if model == "uptake"
    tia = a(1) * (1/rates(1) - 1/rates(2));
else
    tia = sum(a(:) ./ rates(:));
end
end

function w = windowIntegral(model, a, rates, t1, t2)
% Exact integral of the fitted curve over [t1, t2].
seg = @(r) (exp(-r*t1) - exp(-r*t2)) ./ r;
if model == "uptake"
    w = a(1) * (seg(rates(1)) - seg(rates(2)));
else
    w = sum(a(:) .* seg(rates(:)));
end
end

function msgs = qualityWarnings(fit, lamP, warnFrac)
msgs = strings(0, 1);
if fit.exactlyDetermined
    msgs(end+1, 1) = "Model is exactly determined (as many parameters as time points): zero residual, so fit quality cannot be assessed.";
end
if isfinite(fit.extrapolatedFraction) && fit.extrapolatedFraction > warnFrac
    msgs(end+1, 1) = sprintf("%.0f%% of the time-integrated activity lies outside the measured window [%g, %g] h and is extrapolated.", ...
        100*fit.extrapolatedFraction, fit.t(1), fit.t(end));
end
if fit.lambda(1) - lamP <= 1.001e-6
    msgs(end+1, 1) = "No biological clearance could be resolved (slowest rate sits at the physical-decay bound); TIA is limited only by physical decay.";
end
if numel(fit.lambda) > 1 && any(diff(fit.lambda) >= 19.9)
    msgs(end+1, 1) = "A fitted rate is at the upper search limit: that component is faster than the sampling can resolve.";
end
if all(fit.amplitude == 0)
    msgs(end+1, 1) = "All fitted amplitudes are zero.";
end
if fit.weights == "relative" && ~fit.exactlyDetermined && fit.weightedRMSE > 0.25
    msgs(end+1, 1) = sprintf("Poor fit: typical misfit is %.0f%% of the measured value.", 100*fit.weightedRMSE);
end
end

function th = searchRates(model, obj)
if model == "mono"
    th = fminbnd(@(x) obj(x), log(1e-6), log(5), optimset('TolX', 1e-10));
    return
end
opt = optimset('Display', 'off', 'TolX', 1e-7, 'TolFun', 1e-14, ...
               'MaxFunEvals', 800, 'MaxIter', 800);
switch model
    case {"biexp","uptake"}
        starts = startGrid(log([0.002 0.01 0.04]), log([0.05 0.5]));
    case "triexp"
        starts = startGrid(log([0.003 0.02]), log([0.05 0.4]), log([0.5 3]));
end
best = inf;
th   = starts(1, :);
for i = 1:size(starts, 1)
    [x, fv] = fminsearch(obj, starts(i, :), opt);
    if fv < best
        best = fv;
        th   = x;
    end
end
end

function S = startGrid(varargin)
% All combinations of the given start values, one combination per row.
g = cell(1, nargin);
[g{:}] = ndgrid(varargin{:});
S = cell2mat(cellfun(@(x) x(:), g, 'UniformOutput', false));
end
