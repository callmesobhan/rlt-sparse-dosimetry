function D = absorbedDose(tia, massKg, opts)
%ABSORBEDDOSE Absorbed dose from time-integrated activity (local deposition model).
%   D = rlt.absorbedDose(tia, massKg)             % tia in MBq*h, D in Gy
%   D = rlt.absorbedDose(tia, massKg, Isotope="Y90", AbsorbedFraction=0.9)
%   D = rlt.absorbedDose(tia, massKg, SValue=s)   % s in Gy/(MBq*h)
%
%   Default model: all beta energy is absorbed in the source region,
%       D = tia * decays/(MBq*h) * E_beta * 1.602e-13 J/MeV * phi / mass
%   with phi = AbsorbedFraction (default 1). Photon and cross-organ
%   contributions are ignored; for small lesions, where beta escape matters,
%   pass phi < 1 or a sphere S-value via SValue.
%
%   This is a screening model for comparing sampling schedules and fits. It is
%   not a replacement for MIRD S-values or Monte Carlo dosimetry.
arguments
    tia double
    massKg double {mustBePositive}
    opts.Isotope (1,1) string = "Lu177"
    opts.AbsorbedFraction (1,1) double = 1
    opts.SValue double = []
end
if opts.AbsorbedFraction < 0 || opts.AbsorbedFraction > 1
    error('rlt:badFraction', 'AbsorbedFraction must be in [0, 1].');
end

if isempty(opts.SValue)
    iso          = rlt.isotope(opts.Isotope);
    MeV_to_J     = 1.602176634e-13;
    decaysPerMBqh = 1e6 * 3600;
    S = decaysPerMBqh * iso.meanBetaMeV * MeV_to_J * opts.AbsorbedFraction ./ massKg;
else
    S = opts.SValue;
end
D = tia .* S;
end
