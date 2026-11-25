function res = predictCycle2(tSparse, ASparse, prior, opts)
%PREDICTCYCLE2 Empirical Bayes MAP TIA estimate from 1 or 2 sparse Cycle-2+ measurements.
%   res = rlt.nlme.predictCycle2(tSparse, ASparse, prior)
%   res = rlt.nlme.predictCycle2(tSparse, ASparse, prior, Cycle1Fit=c1Fit)
%
%   Combines sparse Cycle-2+ data (even a single scan or blood draw) with:
%     1. The clinic's empirical Cycle-1 population prior `prior` (from rlt.nlme.fitPopPrior), and
%     2. Optionally, the same patient's Cycle-1 multi-time-point fit `Cycle1Fit`, broadened by
%        physiological inter-cycle variability (`InterCycleAmpCV`, `InterCycleRateCV`).
%
%   Solves the exact log-parameter MAP objective with physical decay constraint
%   lambda_eff = lambda_phys + exp(theta_2), computes the Laplace posterior covariance
%   via central finite-difference Hessian, and returns posterior TIA, 95% credible
%   interval, and comparison against naive Hanscheid single-time-point dosimetry.
%
%   Base MATLAB only; no toolboxes required.
arguments
    tSparse double {mustBeVector, mustBeNonnegative, mustBeFinite}
    ASparse double {mustBeVector, mustBePositive, mustBeFinite}
    prior (1,1) struct
    opts.Cycle1Fit = []
    opts.MeasurementCV (1,1) double = NaN
    opts.InterCycleAmpCV (1,1) double {mustBePositive} = 0.25
    opts.InterCycleRateCV (1,1) double {mustBePositive} = 0.15
    opts.DosePerTIA (1,1) double = NaN
    opts.DoseLimitGy (1,1) double = NaN
end

t = tSparse(:);
A = ASparse(:);
if numel(t) ~= numel(A) || isempty(t)
    error('rlt:badInput', 'tSparse and ASparse must be non-empty vectors of equal length.');
end

lamP = prior.lambdaPhys;
if isnan(opts.MeasurementCV)
    sigObs = prior.sigmaObs;
else
    if opts.MeasurementCV <= 0
        error('rlt:badCV', 'MeasurementCV must be positive.');
    end
    sigObs = sqrt(log(1 + opts.MeasurementCV^2));
end

% Build effective prior (population prior, fused with patient's Cycle-1 posterior if given)
invOmegaPop = inv(prior.omega);
precPrior = invOmegaPop;
rhsPrior  = invOmegaPop * prior.mu;
usedCycle1 = false;

if ~isempty(opts.Cycle1Fit)
    usedCycle1 = true;
    f1 = opts.Cycle1Fit;
    th1 = [log(max(f1.amplitude(1), 1e-12)); log(max(f1.lambda(1) - lamP, 1e-5))];
    sigA = sqrt(log(1 + opts.InterCycleAmpCV^2));
    sigL = sqrt(log(1 + opts.InterCycleRateCV^2));
    omegaC1 = diag([sigA^2, sigL^2]);
    invOmegaC1 = inv(omegaC1);
    precPrior = precPrior + invOmegaC1;
    rhsPrior  = rhsPrior + invOmegaC1 * th1;
end
covPrior = inv(precPrior);
muPrior  = covPrior * rhsPrior;

% Negative log-posterior in log-space theta = [log(A0); log(lambda_bio)]
logObs = log(A);
negLogPost = @(th) mapObjective(th, t, logObs, lamP, sigObs, muPrior, precPrior);

% Multi-start Nelder-Mead around prior mean and observation-implied amplitude
lam0 = lamP + exp(muPrior(2));
a0Guess = mean(A .* exp(lam0 * t));
starts = [
    muPrior';
    log(max(a0Guess, 1e-12)), muPrior(2);
    log(max(a0Guess, 1e-12)), muPrior(2) - 0.5;
    log(max(a0Guess, 1e-12)), muPrior(2) + 0.5
];
optNM = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-10, 'MaxIter', 800, 'MaxFunEvals', 1200);
bestVal = inf;
thMap = muPrior;
for k = 1:size(starts, 1)
    [x, fv] = fminsearch(negLogPost, starts(k, :)', optNM);
    if fv < bestVal
        bestVal = fv;
        thMap = x(:);
    end
end

% Laplace posterior covariance via numerical Hessian at MAP
H = numericalHessian(negLogPost, thMap);
[V, D] = eig(0.5 * (H + H'));
dVals = max(diag(D), 1e-6);
covPost = V * diag(1 ./ dVals) * V';

aMap   = exp(thMap(1));
lamBio = exp(thMap(2));
lamEff = lamP + lamBio;
tiaMap = aMap / lamEff;

% Delta-method variance of log(TIA) = th1 - log(lamP + exp(th2))
gradLogTIA = [1; -lamBio / lamEff];
varLogTIA  = max(gradLogTIA' * covPost * gradLogTIA, 1e-12);
sdLogTIA   = sqrt(varLogTIA);
ci95       = tiaMap * exp([-1.9599639845, 1.9599639845] * sdLogTIA);

% Compare with naive single-time-point Hanscheid estimate (using latest time point)
[~, idxLatest] = max(t);
tiaNaiveSTP = rlt.tiaSinglePoint(t(idxLatest), A(idxLatest), prior.lambdaPop);

res.tia             = tiaMap;
res.tiaCI95         = ci95;
res.tiaRelSD        = sqrt(exp(varLogTIA) - 1);
res.amplitude       = aMap;
res.lambdaEff       = lamEff;
res.lambdaBio       = lamBio;
res.thetaMAP        = thMap;
res.posteriorCov    = covPost;
res.usedCycle1Prior = usedCycle1;
res.tiaNaiveSTP     = tiaNaiveSTP;
res.doseGy          = NaN;
res.doseCI95        = [NaN, NaN];
res.probExceedLimit = NaN;
if ~isnan(opts.DosePerTIA)
    res.doseGy   = tiaMap * opts.DosePerTIA;
    res.doseCI95 = ci95 * opts.DosePerTIA;
    if ~isnan(opts.DoseLimitGy) && opts.DoseLimitGy > 0
        % Log-normal tail probability P(Dose >= DoseLimitGy)
        z = (log(opts.DoseLimitGy / max(res.doseGy, 1e-15))) / sdLogTIA;
        res.probExceedLimit = 0.5 * erfc(z / sqrt(2));
    end
end
res.provenance = rlt.provenance();
end

% ---------------------------------------------------------------------------

function val = mapObjective(th, t, logObs, lamP, sigObs, muPrior, precPrior)
th1 = min(max(th(1), -40), 40);
th2 = min(max(th(2), log(1e-6)), log(20));
lamEff = lamP + exp(th2);
logPred = th1 - lamEff * t;
ll = 0.5 * sum(((logObs - logPred) / sigObs).^2);
diff = [th1; th2] - muPrior;
lp = 0.5 * (diff' * precPrior * diff);
val = ll + lp;
end

function H = numericalHessian(fn, x)
h = 1e-4;
n = numel(x);
H = zeros(n, n);
f0 = fn(x);
for i = 1:n
    ei = zeros(n, 1); ei(i) = h;
    fp = fn(x + ei);
    fm = fn(x - ei);
    H(i, i) = (fp - 2*f0 + fm) / (h^2);
    for j = (i+1):n
        ej = zeros(n, 1); ej(j) = h;
        fpp = fn(x + ei + ej);
        fpm = fn(x + ei - ej);
        fmp = fn(x - ei + ej);
        fmm = fn(x - ei - ej);
        H(i, j) = (fpp - fpm - fmp + fmm) / (4 * h^2);
        H(j, i) = H(i, j);
    end
end
end
