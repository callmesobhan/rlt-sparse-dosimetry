function out = propagateTIA(t, A, opts)
%PROPAGATETIA Propagate assumed measurement error into the TIA estimate.
%   out = rlt.propagateTIA(t, A, CV=0.10, Method="select")
%
%   Each Monte Carlo draw multiplies every measured activity by independent
%   mean-one lognormal noise with coefficient of variation CV (quantification
%   error from calibration, recovery, segmentation), then re-runs the chosen
%   estimator. The spread of the results is the TIA uncertainty caused by
%   measurement error alone, for THIS sampling schedule. It does not include
%   model-form error; for that, compare methods (rlt.estimate) or run
%   rlt.evaluateSchedules on a population.
%
%   Name-value options
%     Method     any method of rlt.estimate (default "select")
%     Isotope    default "Lu177"
%     CV         per-point relative error, default 0.10
%     N          number of draws, default 500
%     Seed       RNG seed (a private stream is used; global RNG untouched)
%     LambdaPop  needed for Method="stp"
%
%   Output struct: tia (N draws), central (estimate on the data as given),
%   median, ci95 ([2.5% 97.5%]), relHalfWidth (half the 95% interval over the
%   median), nFailed.
arguments
    t double {mustBeVector, mustBeNonnegative, mustBeFinite}
    A double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Method (1,1) string = "select"
    opts.Isotope (1,1) string = "Lu177"
    opts.CV (1,1) double {mustBePositive} = 0.10
    opts.N (1,1) double {mustBeInteger, mustBePositive} = 500
    opts.Seed (1,1) double = 1
    opts.LambdaPop (1,1) double = NaN
end
t = t(:);
A = A(:);
n = numel(t);

central = rlt.estimate(opts.Method, t, A, Isotope=opts.Isotope, LambdaPop=opts.LambdaPop);

sigma  = sqrt(log(1 + opts.CV^2));
stream = RandStream('mt19937ar', 'Seed', opts.Seed);
noise  = exp(sigma * randn(stream, opts.N, n) - sigma^2/2);

draws = nan(opts.N, 1);
for i = 1:opts.N
    try
        draws(i) = rlt.estimate(opts.Method, t, A .* noise(i, :)', ...
            Isotope=opts.Isotope, LambdaPop=opts.LambdaPop);
    catch
        draws(i) = NaN;
    end
end

q = rlt.internal.quantileLinear(draws, [0.025 0.5 0.975]);
out.tia          = draws;
out.central      = central;
out.median       = q(2);
out.ci95         = [q(1) q(3)];
out.relHalfWidth = (q(3) - q(1)) / (2*q(2));
out.nFailed      = sum(isnan(draws));
end
