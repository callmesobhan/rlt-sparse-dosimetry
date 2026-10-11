function ts = optimalSinglePointTime(lambdaEff)
%OPTIMALSINGLEPOINTTIME Scan time that minimises sensitivity to the assumed rate.
%   ts = rlt.optimalSinglePointTime(lambdaEff) returns 1/lambdaEff (h), i.e.
%   T_eff/ln(2). For a mono-exponential curve this is where
%   d(relative TIA error)/d(lambdaEff) = 0 (see rlt.tiaSinglePoint).
arguments
    lambdaEff double {mustBePositive}
end
ts = 1 ./ lambdaEff;
end
