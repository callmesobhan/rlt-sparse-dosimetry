function prior = fitPopPrior(series, opts)
%FITPOPPRIOR Fit an empirical log-normal population-PK prior from Cycle-1 patients.
%   prior = rlt.nlme.fitPopPrior(series)
%   prior = rlt.nlme.fitPopPrior(series, Isotope="Lu177", MinPatients=3)
%
%   Commercial RLT dosimetry software runs single-patient fits or fixed literature
%   single-time-point rates. This function ingests a clinic's own Cycle-1
%   multi-time-point cohort and estimates an empirical hierarchical log-normal
%   population prior over [log(A0), log(lambda_bio)] (mono-exponential effective
%   clearance constrained by physical decay lambda_eff = lambda_phys + lambda_bio).
%
%   Inputs
%     series : cell array of structs, each with fields .t (hours) and .A (activity
%              or concentration, NOT decay-corrected).
%
%   Name-value options
%     Isotope      "Lu177" (default), "Y90", "I131"
%     MinPatients  minimum valid Cycle-1 fits required (default 3)
%     CovFloor     minimum variance floor on log-parameters for numerical stability
%                  in small cohorts (default 1e-3)
%
%   Output struct `prior` fields
%     isotope, lambdaPhys, nPatients, mu (2x1 mean of [log(A0); log(lamBio)]),
%     omega (2x2 between-patient covariance), sigmaObs (estimated log-residual SD),
%     lambdaPop (geometric-mean effective clearance rate, 1/h),
%     patientParams (N x 2 matrix of fitted [log(A0), log(lamBio)]),
%     patientTIA (N x 1 vector of Cycle-1 TIAs), fits (cell array of fitTAC structs),
%     provenance.
%
%   Base MATLAB only; no toolboxes required.
arguments
    series cell
    opts.Isotope (1,1) string = "Lu177"
    opts.MinPatients (1,1) double {mustBeInteger, mustBePositive} = 3
    opts.CovFloor (1,1) double {mustBeNonnegative} = 1e-3
end

iso = rlt.isotope(opts.Isotope);
lamP = iso.lambdaPhys;
nIn = numel(series);
if nIn < opts.MinPatients
    error('rlt:tooFewPatients', ...
        'Need at least %d Cycle-1 patient series to fit a population prior (got %d).', ...
        opts.MinPatients, nIn);
end

logA = nan(nIn, 1);
logLamBio = nan(nIn, 1);
tias = nan(nIn, 1);
fits = cell(nIn, 1);
resList = [];

for i = 1:nIn
    s = series{i};
    if ~isstruct(s) || ~isfield(s, 't') || ~isfield(s, 'A')
        error('rlt:badSeries', 'series{%d} must be a struct with fields .t and .A.', i);
    end
    f = rlt.fitTAC(s.t, s.A, Isotope=opts.Isotope, Model="mono");
    fits{i} = f;
    a0 = max(f.amplitude(1), 1e-12);
    lamBio = max(f.lambda(1) - lamP, 1e-5);
    logA(i) = log(a0);
    logLamBio(i) = log(lamBio);
    tias(i) = f.tia;

    pred = max(f.predict(s.t), 1e-15);
    obs = max(s.A(:), 1e-15);
    if numel(obs) > 2
        resList = [resList; log(obs ./ pred)]; %#ok<AGROW>
    end
end

theta = [logA, logLamBio];
mu = mean(theta, 1)';
if nIn > 1
    omega = cov(theta);
else
    omega = zeros(2, 2);
end
% Guarantee strictly positive-definite covariance in tiny cohorts
omega = omega + opts.CovFloor * eye(2);

if numel(resList) >= 3
    sigmaObs = max(std(resList), 0.03);
else
    sigmaObs = 0.10;
end

prior.isotope       = iso.name;
prior.lambdaPhys    = lamP;
prior.nPatients     = nIn;
prior.mu            = mu;
prior.omega         = omega;
prior.sigmaObs      = sigmaObs;
prior.lambdaPop     = lamP + exp(mu(2));
prior.patientParams = theta;
prior.patientTIA    = tias;
prior.fits          = fits;
prior.provenance    = rlt.provenance();
end
