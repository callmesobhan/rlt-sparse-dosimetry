function S = scaleSValue(S0, mRef, mPatient, opts)
%SCALESVALUE Scale a phantom S-value to a patient mass.
%   S = rlt.marrow.scaleSValue(S0, mRef, mPatient)
%   S = rlt.marrow.scaleSValue(S0, mRef, mPatient, Exponent=1)
%
%   S = S0 * (mRef/mPatient)^Exponent. Exponent 1 is the standard first-order
%   rule for self-irradiation by non-penetrating radiation (S proportional to
%   1/mass). For photon cross-dose from the remainder of the body it is only an
%   approximation, commonly applied with total-body mass (see Wessels et al.,
%   J Nucl Med 2004, for a standardised treatment). Treat the choice as a
%   documented assumption, and prefer S-values computed for the patient's size.
arguments
    S0 double {mustBeNonnegative, mustBeFinite}
    mRef (1,1) double {mustBePositive}
    mPatient (1,1) double {mustBePositive}
    opts.Exponent (1,1) double = 1
end
S = S0 .* (mRef / mPatient) ^ opts.Exponent;
end
