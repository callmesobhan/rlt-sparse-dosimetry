function audit = auditBloodProxy(tBlood, CBlood, tMarrow, CMarrow, opts)
%AUDITBLOODPROXY Detect when blood-based red-marrow dosimetry violates the Sgouros proxy assumption.
%   audit = rlt.marrow.auditBloodProxy(tBlood, CBlood, tMarrow, CMarrow, Hct=0.45)
%
%   The blood-based marrow method assumes no specific marrow binding and no skeletal
%   tumour infiltration (constant RMBLR = RMECFF / (1 - Hct)). In 177Lu-DOTATATE and
%   bone-metastatic 177Lu-PSMA patients, specific marrow uptake or lesion spillover makes
%   blood-based marrow dose a severe underestimate.
%
%   This audit runs three independent physical checks on paired blood and vertebral-VOI curves:
%     1. Dose Discrepancy Ratio: R_D = D_image / D_blood (flags VIOLATED if > MaxDoseRatio).
%     2. Effective Clearance Divergence: marrow retains activity slower than blood
%        (lambda_marrow / lambda_blood < MinRateRatio).
%     3. Rising Marrow-to-Blood Concentration Ratio over time (positive slope of
%        C_marrow(t)/C_blood(t), smoking gun of active receptor binding or marrow infiltration).
%
%   Output struct `audit` fields:
%     status ("VALID" | "CAUTION" | "VIOLATED"), doseBloodGy, doseImageGy, doseRatio,
%     lambdaBlood, lambdaMarrow, rateRatio, rmblrSlopePerDay, flags, recommendation,
%     bloodResult, marrowResult, provenance.
arguments
    tBlood double {mustBeVector, mustBeNonnegative, mustBeFinite}
    CBlood double {mustBeVector, mustBeNonnegative, mustBeFinite}
    tMarrow double {mustBeVector, mustBeNonnegative, mustBeFinite}
    CMarrow double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Hct (1,1) double {mustBePositive} = 0.45
    opts.RMECFF (1,1) double {mustBePositive} = 0.19
    opts.Isotope (1,1) string = "Lu177"
    opts.ConcentrationUnit (1,1) string {mustBeMember(opts.ConcentrationUnit, ["MBq/mL","MBq/g"])} = "MBq/g"
    opts.MarrowDensityGPerML (1,1) double = 1.03
    opts.RecoveryCoefficient (1,1) double {mustBePositive} = 1.0
    opts.MaxDoseRatio (1,1) double {mustBePositive} = 1.35
    opts.MinRateRatio (1,1) double {mustBePositive} = 0.75
end

bloodRes = rlt.marrow.selfDose(tBlood, CBlood, Source="blood", ...
    ConcentrationUnit=opts.ConcentrationUnit, Hct=opts.Hct, RMECFF=opts.RMECFF, ...
    Isotope=opts.Isotope);
if opts.ConcentrationUnit == "MBq/mL"
    marrowRes = rlt.marrow.selfDose(tMarrow, CMarrow, Source="marrow", ...
        ConcentrationUnit="MBq/mL", DensityGPerML=opts.MarrowDensityGPerML, ...
        RecoveryCoefficient=opts.RecoveryCoefficient, LesionFree=true, Isotope=opts.Isotope);
else
    marrowRes = rlt.marrow.selfDose(tMarrow, CMarrow, Source="marrow", ...
        ConcentrationUnit="MBq/g", RecoveryCoefficient=opts.RecoveryCoefficient, ...
        LesionFree=true, Isotope=opts.Isotope);
end

doseRatio = marrowRes.doseGy / max(bloodRes.doseGy, 1e-15);

% Terminal effective clearance rates (slowest component)
lamBlood  = bloodRes.fit.lambda(1);
lamMarrow = marrowRes.fit.lambda(1);
rateRatio = lamMarrow / max(lamBlood, 1e-12);

% Evaluate fitted marrow-to-blood concentration ratio across overlapping time window
tLo = max(min(tBlood), min(tMarrow));
tHi = min(max(tBlood), max(tMarrow));
if tHi > tLo
    tGrid = linspace(tLo, tHi, 20)';
    cB = max(bloodRes.fit.predict(tGrid), 1e-15);
    cM = max(marrowRes.fit.predict(tGrid), 1e-15) / opts.RecoveryCoefficient;
    ratioCurve = cM ./ cB;
    % Linear slope of relative ratio per 24 h
    p = polyfit(tGrid / 24, ratioCurve / max(ratioCurve(1), 1e-12), 1);
    rmblrSlopePerDay = p(1);
else
    rmblrSlopePerDay = NaN;
end

flags = strings(0, 1);
nViol = 0;
if doseRatio > opts.MaxDoseRatio
    nViol = nViol + 1;
    flags(end+1, 1) = sprintf("Image-based marrow dose is %.2fx blood-based dose (exceeds %.2fx threshold): blood proxy underestimates marrow dose.", ...
        doseRatio, opts.MaxDoseRatio);
elseif doseRatio < 1 / opts.MaxDoseRatio
    flags(end+1, 1) = sprintf("Image-based marrow dose is %.2fx blood-based dose (unusually low; check VOI partial-volume recovery coefficient).", ...
        doseRatio);
end

if rateRatio < opts.MinRateRatio
    nViol = nViol + 1;
    flags(end+1, 1) = sprintf("Marrow terminal clearance (%.4f /h) is %.0f%% slower than blood (%.4f /h): indicates specific marrow retention.", ...
        lamMarrow, 100*(1 - rateRatio), lamBlood);
end

if isfinite(rmblrSlopePerDay) && rmblrSlopePerDay > 0.15
    nViol = nViol + 1;
    flags(end+1, 1) = sprintf("Marrow-to-blood concentration ratio rises by +%.0f%%/day: violates constant-RMBLR assumption.", ...
        100 * rmblrSlopePerDay);
end

if nViol >= 2 || doseRatio > 1.6
    status = "VIOLATED";
    rec = "Blood-based marrow proxy is INVALID for this patient (specific marrow uptake or infiltration); use vertebral-VOI image-based dosimetry for Cycle 2+.";
elseif nViol == 1 || ~isempty(flags)
    status = "CAUTION";
    rec = "Borderline blood-vs-image marrow discrepancy; verify lesion-free VOI placement and retain an image-based check in Cycle 2.";
else
    status = "VALID";
    rec = "Blood and marrow-VOI kinetics agree with constant-RMBLR assumptions; blood-based sparse monitoring is supported.";
end

audit.status           = status;
audit.doseBloodGy      = bloodRes.doseGy;
audit.doseImageGy      = marrowRes.doseGy;
audit.doseRatio        = doseRatio;
audit.lambdaBlood      = lamBlood;
audit.lambdaMarrow     = lamMarrow;
audit.rateRatio        = rateRatio;
audit.rmblrSlopePerDay = rmblrSlopePerDay;
audit.flags            = flags;
audit.recommendation   = rec;
audit.bloodResult      = bloodRes;
audit.marrowResult     = marrowRes;
audit.provenance       = rlt.provenance();
end
