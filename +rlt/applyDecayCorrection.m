function Acorr = applyDecayCorrection(t, A, isotope)
%APPLYDECAYCORRECTION Decay-correct measured activities to the reference time (t = 0).
%   Acorr = rlt.applyDecayCorrection(t, A, "Lu177"), the inverse of
%   rlt.undoDecayCorrection. Use it only to present data; fit uncorrected data.
arguments
    t double {mustBeNonnegative, mustBeFinite}
    A double {mustBeNonnegative, mustBeFinite}
    isotope (1,1) string = "Lu177"
end
if ~isequal(size(t), size(A))
    error('rlt:sizeMismatch', 't and A must have the same size.');
end
Acorr = A .* exp(rlt.isotope(isotope).lambdaPhys .* t);
end
