function res = totalDose(contrib)
%TOTALDOSE Sum red-marrow dose contributions from several source regions.
%   res = rlt.marrow.totalDose(contrib)
%
%   contrib is a table with one row per source region and columns
%     Source        text label, unique ("Marrow self", "Kidneys", "RemainderOfBody", ...)
%     TIA_MBqh      time-integrated activity of that source (MBq h)
%     S_GyPerMBqh   S-value from that source to the red marrow (Gy per MBq h)
%   Total marrow dose = sum over sources of TIA x S. A typical set: marrow self
%   (blood or VOI based), remainder of body, each major organ, each tumour
%   (Gosewisch et al., EJNMMI Res 2018). Every term must be listed explicitly;
%   rlt.io.lookupS errors on a missing S-value so nothing silently counts as zero.
%
%   res.table   per-source TIA, S, dose (Gy) and share of the total (%)
%   res.totalGy total marrow dose in Gy
arguments
    contrib table
end
need = ["Source","TIA_MBqh","S_GyPerMBqh"];
if ~all(ismember(need, string(contrib.Properties.VariableNames)))
    error('rlt:badContrib', 'contrib needs columns: %s.', strjoin(need, ", "));
end
if height(contrib) == 0
    error('rlt:badContrib', 'contrib has no rows.');
end
src = string(contrib.Source);
tia = double(contrib.TIA_MBqh);
sv  = double(contrib.S_GyPerMBqh);
if numel(unique(src)) ~= numel(src)
    error('rlt:badContrib', 'Source labels must be unique (a source counted twice would double the dose).');
end
if any(~isfinite(tia) | tia < 0) || any(~isfinite(sv) | sv < 0)
    error('rlt:badContrib', 'TIA_MBqh and S_GyPerMBqh must be finite and non-negative.');
end

D = tia .* sv;
tot = sum(D);
if tot > 0
    share = 100 * D / tot;
else
    share = zeros(size(D));
end
res.table = table(src, tia, sv, D, share, ...
    'VariableNames', {'Source','TIA_MBqh','S_GyPerMBqh','DoseGy','SharePct'});
res.totalGy = tot;
end
