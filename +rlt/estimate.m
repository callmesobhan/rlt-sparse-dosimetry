function tia = estimate(method, t, A, opts)
%ESTIMATE Single entry point for every TIA estimator in the toolbox.
%   tia = rlt.estimate(method, t, A, Isotope=..., LambdaPop=..., Origin=..., Models=...)
%
%   method
%     "mono" | "biexp" | "triexp" | "uptake"   fit that model (rlt.fitTAC)
%     "select"                    AICc selection among Models (rlt.selectModel)
%     "hybrid"                    trapezoid + tail (rlt.tiaHybrid, uses Origin)
%     "stp"                       single time point (rlt.tiaSinglePoint). Uses
%                                 the scan closest to 1/LambdaPop; LambdaPop
%                                 (1/h) is required.
%
%   Origin  passed to rlt.tiaHybrid (default "zero"; use "decay" for blood)
%   Models  candidate models for "select" (default mono, biexp, uptake)
arguments
    method (1,1) string {mustBeMember(method, ["mono","biexp","triexp","uptake","select","hybrid","stp"])}
    t double {mustBeVector}
    A double {mustBeVector}
    opts.Isotope (1,1) string = "Lu177"
    opts.LambdaPop (1,1) double = NaN
    opts.Origin (1,1) string = "zero"
    opts.Models (1,:) string = ["mono","biexp","uptake"]
end
t = t(:);
A = A(:);

switch method
    case {"mono","biexp","triexp","uptake"}
        tia = rlt.fitTAC(t, A, Isotope=opts.Isotope, Model=method).tia;
    case "select"
        tia = rlt.selectModel(t, A, Isotope=opts.Isotope, Models=opts.Models).tia;
    case "hybrid"
        tia = rlt.tiaHybrid(t, A, Isotope=opts.Isotope, Origin=opts.Origin);
    case "stp"
        if isnan(opts.LambdaPop)
            error('rlt:missingLambdaPop', 'Method "stp" needs LambdaPop (population effective rate, 1/h).');
        end
        lamEff = max(opts.LambdaPop, rlt.isotope(opts.Isotope).lambdaPhys);
        [~, k] = min(abs(t - rlt.optimalSinglePointTime(lamEff)));
        tia = rlt.tiaSinglePoint(t(k), A(k), lamEff);
end
end
