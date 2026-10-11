function checkSeries(t, A)
%CHECKSERIES Validate a time-activity series (column vectors, times in h).
if numel(t) ~= numel(A)
    error('rlt:badSeries', 't and A must have the same number of elements.');
end
if numel(t) < 2
    error('rlt:tooFewPoints', 'At least 2 time points are required.');
end
if any(diff(t) <= 0)
    error('rlt:badSeries', 'Time points must be strictly increasing.');
end
if ~any(A > 0)
    error('rlt:badSeries', 'All activities are zero.');
end
end
