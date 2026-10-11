function A = undoDecayCorrection(t, Acorr, isotope)
%UNDODECAYCORRECTION Turn decay-corrected activities back into measured activities.
%   A = rlt.undoDecayCorrection(t, Acorr, "Lu177")
%
%   Many quantification tools report activity decay-corrected to the injection
%   (or scan-start) time: Acorr(t) = A(t)*exp(+lambda_phys*t). rlt.fitTAC models
%   physical decay inside its rates, so it must be given A(t) = Acorr(t)*exp(-lambda_phys*t).
%   Feeding decay-corrected data to the fits silently overestimates the TIA.
%   t must be hours after the reference time of the correction.
arguments
    t double {mustBeNonnegative, mustBeFinite}
    Acorr double {mustBeNonnegative, mustBeFinite}
    isotope (1,1) string = "Lu177"
end
if ~isequal(size(t), size(Acorr))
    error('rlt:sizeMismatch', 't and Acorr must have the same size.');
end
A = Acorr .* exp(-rlt.isotope(isotope).lambdaPhys .* t);
end
