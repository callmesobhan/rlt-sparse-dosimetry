function e = elasticity(t, C, parameter, base, varargin)
%ELASTICITY Sensitivity of the blood-based self-dose to one input: d ln(D) / d ln(p).
%   e = rlt.marrow.elasticity(t, C, "Hct", 0.45)
%   e = rlt.marrow.elasticity(t, C, "RMECFF", 0.19, "Hct", 0.45)
%   e = rlt.marrow.elasticity(t, C, "AbsorbedFraction", 0.8, "Hct", 0.45)
%
%   Central finite difference of rlt.marrow.selfDose with the parameter moved by
%   +/-5% around base. An elasticity of 1 means a 10% error in the parameter gives
%   a 10% error in dose. Analytic values for the blood basis: RMECFF -> 1,
%   Hct -> Hct/(1-Hct), AbsorbedFraction -> 1. Extra name-value pairs (as
%   "Name", value) are passed to selfDose. parameter is "Hct", "RMECFF" or
%   "AbsorbedFraction"; for AbsorbedFraction keep base <= 0.95 so base*1.05 stays <= 1.
parameter = string(parameter);
if ~ismember(parameter, ["Hct","RMECFF","AbsorbedFraction"])
    error('rlt:badParameter', 'parameter must be "Hct", "RMECFF" or "AbsorbedFraction".');
end
if ~isscalar(base) || ~isfinite(base) || base <= 0
    error('rlt:badParameter', 'base must be a positive scalar.');
end
h  = 0.05;
lo = rlt.marrow.selfDose(t, C, varargin{:}, char(parameter), base*(1 - h)).doseGy;
hi = rlt.marrow.selfDose(t, C, varargin{:}, char(parameter), base*(1 + h)).doseGy;
e  = (log(hi) - log(lo)) / (log(1 + h) - log(1 - h));
end
