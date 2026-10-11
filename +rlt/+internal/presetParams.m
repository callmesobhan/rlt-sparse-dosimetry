function P = presetParams(preset)
%PRESETPARAMS Distributions of virtual-patient kinetics.
%
%   *** ILLUSTRATIVE VALUES ONLY. They are NOT measured data. ***
%   Replace them with distributions estimated from your own cohort before
%   drawing any conclusion about a real sampling protocol.
%
%   Fields
%     nComp     number of washout exponentials K
%     frac      1 x K mean fraction of the amplitude in each component, ordered
%               slowest -> fastest
%     fracGsd   geometric SD of the multiplicative scatter on those fractions
%               (fractions are renormalised to sum to 1)
%     Teff      K x 2 [median, geometric SD] of each component's effective
%               half-life (h); the drawn rates are sorted so component 1 is slowest
%     amp       struct(median, gsd): total amplitude at t = 0 (washout phases)
%     hasUptake true if an uptake term is subtracted so that A(0) = 0
%     Tup       struct(median, gsd): uptake half-time (h), used if hasUptake
%     units     what A(t) is measured in
switch string(preset)
    case "kidney"
        P.nComp = 2;
        P.frac = [0.65 0.35];       P.fracGsd = 1.25;
        P.Teff = [45 1.3; 8 1.5];
        P.amp = struct('median', 0.020, 'gsd', 1.4);
        P.hasUptake = true;
        P.Tup = struct('median', 1.0, 'gsd', 1.5);
        P.units = "fraction of injected activity in the organ";
    case "tumour"
        P.nComp = 2;
        P.frac = [0.85 0.15];       P.fracGsd = 1.3;
        P.Teff = [80 1.3; 15 1.5];
        P.amp = struct('median', 0.010, 'gsd', 2.0);
        P.hasUptake = true;
        P.Tup = struct('median', 2.0, 'gsd', 1.5);
        P.units = "fraction of injected activity in the lesion";
    case "blood"
        % Falls from injection; three phases so that mono- and bi-exponential
        % fits show model-form bias against the truth. Amplitude is the
        % concentration at t = 0 as a fraction of injected activity per mL
        % (about 1/blood volume).
        P.nComp = 3;
        P.frac = [0.03 0.17 0.80];  P.fracGsd = 1.3;
        P.Teff = [48 1.4; 8 1.4; 1.2 1.5];
        P.amp = struct('median', 2.0e-4, 'gsd', 1.25);
        P.hasUptake = false;
        P.Tup = struct('median', NaN, 'gsd', NaN);
        P.units = "fraction of injected activity per mL of blood";
    otherwise
        error('rlt:unknownPreset', 'Unknown preset "%s". Use "kidney", "tumour" or "blood".', string(preset));
end
end
