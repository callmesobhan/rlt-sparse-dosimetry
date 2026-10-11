function tia = tiaSinglePoint(ts, As, lambdaEff)
%TIASINGLEPOINT Time-integrated activity from one measurement and a population rate.
%   tia = rlt.tiaSinglePoint(ts, As, lambdaEff)
%
%   Assumes mono-exponential decay with the population effective rate
%   lambdaEff (1/h), so TIA = As * exp(lambdaEff*ts) / lambdaEff.
%   If the true rate is lambdaEff + d, the relative error is zero to first
%   order in d when ts = 1/lambdaEff; see rlt.optimalSinglePointTime.
%   This is the idea behind single-time-point dosimetry for 177Lu (Hanscheid
%   et al., J Nucl Med 2018).
arguments
    ts double {mustBeNonnegative}
    As double {mustBeNonnegative}
    lambdaEff double {mustBePositive}
end
tia = As .* exp(lambdaEff .* ts) ./ lambdaEff;
end
