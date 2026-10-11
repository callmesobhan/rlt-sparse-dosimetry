function q = quantileLinear(x, p)
%QUANTILELINEAR Empirical quantiles with linear interpolation (no toolbox).
%   q = quantileLinear(x, p) returns the quantiles of x at probabilities p
%   (0..1). NaNs in x are ignored. q has the same size as p.
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0
    q = nan(size(p));
    return
end
pos = 1 + (n - 1) * p(:);
lo  = floor(pos);
hi  = ceil(pos);
w   = pos - lo;
q   = x(lo) .* (1 - w) + x(hi) .* w;
q   = reshape(q, size(p));
end
