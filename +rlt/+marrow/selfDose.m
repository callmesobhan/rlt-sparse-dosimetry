function res = selfDose(t, C, opts)
%SELFDOSE Red-marrow self-dose from a blood or marrow-VOI concentration curve.
%   res = rlt.marrow.selfDose(t, C, Hct=0.45)                       % blood-based
%   res = rlt.marrow.selfDose(t, C, Source="marrow", LesionFree=true) % image-based VOI
%
%   t : sample times (h after administration)
%   C : activity concentration, NOT decay-corrected (MBq/mL or MBq/g)
%
%   Blood-based method (AAPM/Sgouros; EANM guideline): the red-marrow activity
%   concentration is a fixed multiple of the blood concentration,
%       C_RM(t) = RMBLR * C_blood(t),   RMBLR = RMECFF/(1 - Hct),
%   fitted with a decreasing exponential model. Under the default local-deposition
%   model the self-dose is
%       D = RMBLR * TIAconc_blood[MBq h/g] * 1000 g/kg * 3.6e9 * E_beta * 1.602e-13 * phi
%   and the marrow mass cancels out.
%
%   WHEN THIS IS NOT VALID (see docs/MARROW.md): the ratio assumes no specific
%   marrow uptake and no marrow involvement. Specific marrow uptake (reported
%   for 177Lu-DOTATATE) and bone metastases make the blood method a lower
%   bound; compare it with the image-based result and with the model spread
%   reported here.
%
%   Name-value options
%     Source           "blood" (default), "plasma", or "marrow" (the curve is already
%                      a marrow/vertebral VOI concentration)
%     ConcentrationUnit "MBq/mL" (default) or "MBq/g"
%     DensityGPerML    tissue density for MBq/mL -> MBq/g. Default 1.06 for blood and
%                      plasma; required for Source="marrow" with MBq/mL.
%     Isotope          default "Lu177"
%     Hct              haematocrit fraction (needed for Source="blood" unless Ratio or RatioFcn is given)
%     RMECFF           default 0.19
%     Ratio            constant marrow-to-source ratio; overrides Hct/RMECFF
%     RatioFcn         handle ratio(t) for a time-dependent ratio (t in h, vectorised);
%                      the literature reports ratios that rise with time
%     Scale, RecoveryCoefficient   for Source="marrow": effective ratio = Scale/RecoveryCoefficient
%                      (Scale: marrow concentration over VOI concentration; default 1.
%                      RecoveryCoefficient: measured/true from phantom data; default 1.)
%     LesionFree       assert that a marrow VOI contains no tumour (default false -> warning)
%     Sex, MarrowMassG marrow mass for reporting the marrow TIA in MBq h
%     Model            "select" (default; AICc among Models) or one model
%     Models           default ["mono","biexp","triexp"]
%     AbsorbedFraction beta fraction absorbed in marrow, default 1 (an upper bound)
%     SelfS, SelfSRefMassG   S(RM<-RM) in Gy/(MBq h) for a marrow of SelfSRefMassG grams;
%                      replaces the local-deposition model (mass-scaled automatically)
%     AdministeredMBq  adds Gy per GBq to the result
%     SpreadWarn       warn if the model-to-model dose spread exceeds this (default 0.25)
%
%   res fields: doseGy, dosePerGBq, marrowTIA (MBq h; NaN without a mass),
%   marrowConcTIA (MBq h per g marrow), ratioEffective, fit (selected model),
%   candidates (every model with its dose), modelSpread, warnings, assumptions,
%   inputs, provenance.
arguments
    t double {mustBeVector, mustBeNonnegative, mustBeFinite}
    C double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Source (1,1) string {mustBeMember(opts.Source, ["blood","plasma","marrow"])} = "blood"
    opts.ConcentrationUnit (1,1) string {mustBeMember(opts.ConcentrationUnit, ["MBq/mL","MBq/g"])} = "MBq/mL"
    opts.DensityGPerML (1,1) double = NaN
    opts.Isotope (1,1) string = "Lu177"
    opts.Hct (1,1) double = NaN
    opts.RMECFF (1,1) double = 0.19
    opts.Ratio double = []
    opts.RatioFcn = []
    opts.Scale (1,1) double {mustBePositive} = 1
    opts.RecoveryCoefficient (1,1) double {mustBePositive} = 1
    opts.LesionFree (1,1) logical = false
    opts.Sex (1,1) string = ""
    opts.MarrowMassG double = []
    opts.Model (1,1) string {mustBeMember(opts.Model, ["select","mono","biexp","triexp"])} = "select"
    opts.Models (1,:) string = ["mono","biexp","triexp"]
    opts.AbsorbedFraction (1,1) double = 1
    opts.SelfS double = []
    opts.SelfSRefMassG double = []
    opts.AdministeredMBq double = []
    opts.SpreadWarn (1,1) double = 0.25
end
t = t(:);
C = C(:);
rlt.internal.checkSeries(t, C);
iso = rlt.isotope(opts.Isotope);
msgs = strings(0, 1);
assumptions = strings(0, 1);

if any(~ismember(opts.Models, ["mono","biexp","triexp"]))
    error('rlt:badModels', 'Models must be a subset of mono, biexp, triexp (decreasing curves).');
end
if opts.AbsorbedFraction < 0 || opts.AbsorbedFraction > 1
    error('rlt:badFraction', 'AbsorbedFraction must be in [0, 1].');
end
if ~isempty(opts.Ratio) && ~isempty(opts.RatioFcn)
    error('rlt:badRatio', 'Give Ratio or RatioFcn, not both.');
end
if ~isempty(opts.Ratio) && ~isscalar(opts.Ratio)
    error('rlt:badRatio', 'Ratio must be a scalar.');
end
if ~isempty(opts.RatioFcn) && ~isa(opts.RatioFcn, 'function_handle')
    error('rlt:badRatio', 'RatioFcn must be a function handle ratio(t).');
end
if ~isempty(opts.SelfS) && (~isscalar(opts.SelfS) || opts.SelfS < 0 || ~isscalar(opts.SelfSRefMassG) && ~isempty(opts.SelfSRefMassG))
    error('rlt:badSelfS', 'SelfS must be a non-negative scalar and SelfSRefMassG a positive scalar.');
end
if ~isempty(opts.SelfS)
    if isempty(opts.SelfSRefMassG)
        error('rlt:needRefMass', 'SelfS needs SelfSRefMassG (the phantom marrow mass behind that S-value).');
    end
    if opts.AbsorbedFraction ~= 1
        msgs(end+1, 1) = "AbsorbedFraction is ignored because SelfS was supplied.";
    end
end

% ---- concentration in MBq per gram ---------------------------------------
if opts.ConcentrationUnit == "MBq/mL"
    dens = opts.DensityGPerML;
    if isnan(dens)
        if opts.Source == "marrow"
            error('rlt:needDensity', 'Source="marrow" with MBq/mL needs DensityGPerML for the VOI tissue.');
        end
        dens = 1.06;
        assumptions(end+1, 1) = "Blood density 1.06 g/mL used to convert MBq/mL to MBq/g.";
    end
    Cg = C / dens;
else
    Cg = C;
end

% ---- marrow-to-source ratio ----------------------------------------------
ratio    = NaN;
ratioFcn = opts.RatioFcn;
switch opts.Source
    case {"blood","plasma"}
        if ~isempty(opts.Ratio)
            ratio = opts.Ratio;
            assumptions(end+1, 1) = sprintf("Marrow-to-%s ratio fixed by the user at %.3f.", opts.Source, ratio);
        elseif ~isempty(ratioFcn)
            assumptions(end+1, 1) = "Time-dependent marrow-to-source ratio supplied by the user (RatioFcn).";
        else
            if opts.Source == "blood"
                if isnan(opts.Hct)
                    error('rlt:needHct', 'Source="blood" needs Hct (or Ratio / RatioFcn).');
                end
                [ratio, m] = rlt.marrow.rmblr(opts.Hct, RMECFF=opts.RMECFF, Basis="blood");
                assumptions(end+1, 1) = sprintf("RMBLR = RMECFF/(1-Hct) = %.2f/(1-%.2f) = %.3f (Sgouros 1993); constant in time.", ...
                    opts.RMECFF, opts.Hct, ratio);
                if strlength(m) > 0, msgs(end+1, 1) = m; end
            else
                ratio = rlt.marrow.rmblr(Basis="plasma", RMECFF=opts.RMECFF);
                assumptions(end+1, 1) = sprintf("Red-marrow-to-plasma ratio = RMECFF = %.3f; constant in time.", ratio);
            end
        end
        assumptions(end+1, 1) = "No specific marrow uptake and no marrow involvement are assumed; otherwise this self-dose is a lower bound.";
    case "marrow"
        if ~isempty(ratioFcn)
            error('rlt:badRatio', 'RatioFcn is not used with Source="marrow".');
        end
        ratio = opts.Scale / opts.RecoveryCoefficient;
        if ~isempty(opts.Ratio), ratio = opts.Ratio; end
        assumptions(end+1, 1) = sprintf("Marrow concentration = VOI concentration x %.3f (Scale/RecoveryCoefficient).", ratio);
        if ~opts.LesionFree
            msgs(end+1, 1) = "VOI not asserted lesion-free (LesionFree=false): tumour or bone-lesion activity inside the VOI can overestimate marrow activity severalfold.";
        end
end
if isempty(ratioFcn)
    if ~isfinite(ratio) || ratio <= 0
        error('rlt:badRatio', 'The marrow-to-source ratio must be a positive finite number.');
    end
else
    rv = ratioFcn(t);
    if ~isequal(size(rv), size(t)) || any(~isfinite(rv)) || any(rv <= 0)
        error('rlt:badRatio', 'RatioFcn(t) must return positive finite values, one per time point.');
    end
end

% ---- marrow mass (reporting only under local deposition) ------------------
massG = NaN;
if ~isempty(opts.MarrowMassG)
    if ~isscalar(opts.MarrowMassG) || ~isfinite(opts.MarrowMassG) || opts.MarrowMassG <= 0
        error('rlt:badMass', 'MarrowMassG must be a positive scalar.');
    end
    massG = opts.MarrowMassG;
    assumptions(end+1, 1) = sprintf("Marrow mass %.0f g supplied by the user (used only to report the marrow TIA).", massG);
elseif strlength(opts.Sex) > 0
    massG = rlt.marrow.referenceMass(opts.Sex);
    assumptions(end+1, 1) = sprintf("Marrow mass %.0f g is the ICRP reference adult value, not patient-specific (used only to report the marrow TIA).", massG);
end

% ---- curve fits -----------------------------------------------------------
if opts.Model == "select"
    [best, ~, fits] = rlt.selectModel(t, Cg, Isotope=opts.Isotope, Models=opts.Models);
else
    fits = {rlt.fitTAC(t, Cg, Isotope=opts.Isotope, Model=opts.Model)};
    best = fits{1};
end
nF = numel(fits);
concTIA = zeros(nF, 1);
doseGy  = zeros(nF, 1);
for i = 1:nF
    concTIA(i) = marrowConcTIA(fits{i}, ratio, ratioFcn, iso);
    doseGy(i)  = doseFromConc(concTIA(i), opts);
end
sel = find(cellfun(@(f) f.model == best.model, fits), 1);

modelName = strings(nF, 1);
nPar = zeros(nF, 1);
tiaSrc = zeros(nF, 1);
aicc = zeros(nF, 1);
for i = 1:nF
    modelName(i) = fits{i}.model;
    nPar(i)  = fits{i}.nParams;
    tiaSrc(i) = fits{i}.tia;
    aicc(i)  = fits{i}.aicc;
end
isSel = false(nF, 1);
isSel(sel) = true;
candidates = table(modelName, nPar, tiaSrc, concTIA, doseGy, aicc, isSel, ...
    'VariableNames', {'Model','nParams','SourceTIA_MBqhPerG','MarrowConcTIA_MBqhPerG','DoseGy','AICc','Selected'});

spread = NaN;
if nF > 1
    spread = (max(doseGy) - min(doseGy)) / median(doseGy);
    if spread > opts.SpreadWarn
        msgs(end+1, 1) = sprintf("Dose differs by %.0f%% between the fitted models (%.3g to %.3g Gy): the blood curve model matters for this patient.", ...
            100*spread, min(doseGy), max(doseGy));
    end
end
for w = best.warnings'
    msgs(end+1, 1) = "Curve fit: " + w;
end
if isempty(opts.SelfS)
    assumptions(end+1, 1) = sprintf("Local deposition: %.0f%% of the mean beta energy (%.3f MeV) is absorbed in the marrow; photons and cross-dose from other regions are NOT included.", ...
        100*opts.AbsorbedFraction, iso.meanBetaMeV);
else
    assumptions(end+1, 1) = sprintf("Self-dose from the supplied S(RM<-RM) = %.3g Gy/(MBq h) for a %.0f g marrow, mass-scaled; other regions are NOT included.", ...
        opts.SelfS, opts.SelfSRefMassG);
end

res.source         = opts.Source;
res.isotope        = iso.name;
res.doseGy         = doseGy(sel);
res.dosePerGBq     = NaN;
if ~isempty(opts.AdministeredMBq)
    if ~isscalar(opts.AdministeredMBq) || ~(opts.AdministeredMBq > 0)
        error('rlt:badActivity', 'AdministeredMBq must be a positive scalar.');
    end
    res.dosePerGBq = res.doseGy / (opts.AdministeredMBq / 1000);
end
res.marrowConcTIA  = concTIA(sel);
res.marrowTIA      = concTIA(sel) * massG;
res.massG          = massG;
res.ratio          = ratio;
res.ratioEffective = concTIA(sel) / best.tia;
res.fit            = best;
res.candidates     = candidates;
res.modelSpread    = spread;
res.warnings       = msgs;
res.assumptions    = assumptions;
res.inputs         = struct('t', t, 'C', C, 'unit', opts.ConcentrationUnit, 'options', opts);
res.provenance     = rlt.provenance();
end

% ---------------------------------------------------------------------------

function c = marrowConcTIA(fit, ratio, ratioFcn, iso)
% Time-integrated marrow activity concentration, MBq h per g of marrow.
if isempty(ratioFcn)
    c = ratio * fit.tia;
else
    f = @(x) reshape(ratioFcn(x(:)), size(x)) .* reshape(fit.predict(x(:)), size(x));
    c = integral(f, 0, 40*iso.halfLifeH, 'RelTol', 1e-8, 'AbsTol', 1e-14*max(fit.A));
end
end

function d = doseFromConc(c, opts)
% Self-dose in Gy from the marrow concentration TIA (MBq h per g marrow).
if isempty(opts.SelfS)
    d = rlt.absorbedDose(c*1000, 1, Isotope=opts.Isotope, AbsorbedFraction=opts.AbsorbedFraction);
else
    % TIA = c*m and S scales as mRef/m, so the patient mass cancels.
    d = c * opts.SelfS * opts.SelfSRefMassG;
end
end
