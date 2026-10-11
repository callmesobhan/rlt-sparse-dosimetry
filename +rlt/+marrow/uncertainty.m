function out = uncertainty(contrib, opts)
%UNCERTAINTY Monte Carlo uncertainty of the total marrow dose.
%   out = rlt.marrow.uncertainty(contrib)
%
%   contrib is the table of rlt.marrow.totalDose plus two columns that YOU must
%   fill, because no default is defensible for every centre:
%     TIA_CV   coefficient of variation of that source's TIA (independent error)
%     S_CV     coefficient of variation of that source's S-value
%   Each draw multiplies TIA and S by independent mean-one lognormal factors.
%
%   Name-value options
%     N          draws (default 5000)
%     Seed       RNG seed (private stream)
%     SharedTIACV  extra CV applied as ONE common factor to every TIA in a draw,
%                  for errors shared by all image-derived terms, such as the
%                  scanner calibration factor (default 0). Ignoring it
%                  understates the uncertainty of the sum.
%
%   out fields
%     total      N x 1 draws of total dose (Gy)
%     median, ci95 ([2.5% 97.5%]), relHalfWidth
%     perSource  table: Source, MedianGy, Lo95Gy, Hi95Gy, VarianceSharePct
%                (approximate: squared correlation of each source's log-dose with
%                the log total, normalised to 100)
arguments
    contrib table
    opts.N (1,1) double {mustBeInteger, mustBePositive} = 5000
    opts.Seed (1,1) double = 1
    opts.SharedTIACV (1,1) double {mustBeNonnegative} = 0
end
need = ["Source","TIA_MBqh","S_GyPerMBqh","TIA_CV","S_CV"];
miss = need(~ismember(need, string(contrib.Properties.VariableNames)));
if ~isempty(miss)
    error('rlt:missingUncertainty', ...
        'contrib is missing column(s): %s. State TIA_CV and S_CV for every source; there are no defaults.', ...
        strjoin(miss, ", "));
end
tia = double(contrib.TIA_MBqh);
sv  = double(contrib.S_GyPerMBqh);
cvT = double(contrib.TIA_CV);
cvS = double(contrib.S_CV);
if any(~isfinite(cvT) | cvT < 0) || any(~isfinite(cvS) | cvS < 0)
    error('rlt:badCV', 'TIA_CV and S_CV must be finite and non-negative.');
end
M = numel(tia);
N = opts.N;
stream = RandStream('mt19937ar', 'Seed', opts.Seed);

lnFactor = @(cv, nRows) exp(sqrt(log(1 + cv.^2)) .* randn(stream, nRows, numel(cv)) ...
                            - log(1 + cv.^2)/2);
fT = lnFactor(cvT(:)', N);                         % N x M
fS = lnFactor(cvS(:)', N);
if opts.SharedTIACV > 0
    fT = fT .* lnFactor(opts.SharedTIACV, N);      % N x 1, common to all sources
end
Dk    = (tia(:)' .* sv(:)') .* fT .* fS;           % N x M
total = sum(Dk, 2);

q = @rlt.internal.quantileLinear;
qt = q(total, [0.025 0.5 0.975]);
out.total        = total;
out.median       = qt(2);
out.ci95         = [qt(1) qt(3)];
out.relHalfWidth = (qt(3) - qt(1)) / (2*qt(2));

share = zeros(M, 1);
lt = log(total);
for k = 1:M
    lk = log(Dk(:, k));
    if all(isfinite(lk)) && std(lk) > 0 && std(lt) > 0
        c = corrcoef(lk, lt);
        share(k) = c(1, 2)^2;
    end
end
if sum(share) > 0
    share = 100 * share / sum(share);
end
medK = zeros(M, 1);  loK = zeros(M, 1);  hiK = zeros(M, 1);
for k = 1:M
    qk = q(Dk(:, k), [0.025 0.5 0.975]);
    loK(k) = qk(1);  medK(k) = qk(2);  hiK(k) = qk(3);
end
out.perSource = table(string(contrib.Source), medK, loK, hiK, share, ...
    'VariableNames', {'Source','MedianGy','Lo95Gy','Hi95Gy','VarianceSharePct'});
end
