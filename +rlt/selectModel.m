function [best, cand, fits] = selectModel(t, A, opts)
%SELECTMODEL Choose among kinetic models by AICc, refusing to over-fit sparse data.
%   [best, cand, fits] = rlt.selectModel(t, A)
%   [...] = rlt.selectModel(t, A, Models=["mono","biexp","triexp"])
%
%   Every requested model with nParams <= numel(t) is fitted. Only models with a
%   finite AICc (numel(t) >= nParams + 2) compete. If none has a finite AICc,
%   which is the usual case with 2-4 scans, the first listed model that was
%   fitted (by default mono-exponential) is returned: the data cannot justify
%   anything more complex.
%
%   best : fit struct (see rlt.fitTAC) of the selected model
%   cand : table of all fitted candidates (Model, nParams, TIA, SSE, AICc, Selected)
%   fits : cell array with every candidate's fit struct, in the order of cand
arguments
    t double {mustBeVector, mustBeNonnegative, mustBeFinite}
    A double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Isotope (1,1) string = "Lu177"
    opts.Models (1,:) string = ["mono","biexp","uptake"]
    opts.Weights (1,1) string {mustBeMember(opts.Weights, ["relative","none"])} = "relative"
end
t = t(:);
A = A(:);
n = numel(t);

fits = {};
for m = opts.Models
    if n >= rlt.internal.nParams(m)
        fits{end+1} = rlt.fitTAC(t, A, Isotope=opts.Isotope, Model=m, Weights=opts.Weights); %#ok<AGROW>
    end
end
if isempty(fits)
    error('rlt:tooFewPoints', 'No requested model can be fitted with %d time points.', n);
end

aicc = cellfun(@(f) f.aicc, fits);
if any(isfinite(aicc))
    aicc(~isfinite(aicc)) = inf;
    [~, k] = min(aicc);
else
    k = 1;
end
best = fits{k};

model = strings(numel(fits), 1);
for i = 1:numel(fits)
    model(i) = fits{i}.model;
end
nParams  = cellfun(@(f) f.nParams, fits)';
tia      = cellfun(@(f) f.tia, fits)';
sse      = cellfun(@(f) f.sse, fits)';
aiccCol  = cellfun(@(f) f.aicc, fits)';
selected = false(numel(fits), 1);
selected(k) = true;
cand = table(model, nParams, tia, sse, aiccCol, selected, ...
    'VariableNames', {'Model','nParams','TIA','SSE','AICc','Selected'});
end
