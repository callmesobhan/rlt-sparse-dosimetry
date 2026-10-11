function m = referenceMass(sex)
%REFERENCEMASS Reference adult red (active) marrow mass in grams.
%   m = rlt.marrow.referenceMass("M")   % 1170 g
%   m = rlt.marrow.referenceMass("F")   % 900 g
%
%   ICRP Publications 70 and 89 reference adult values (as cited in Pichardo
%   et al., J Nucl Med 2007). These are population references, not patient
%   values: real marrow mass varies widely with age, sex, prior therapy and
%   disease. For the self-dose from a local-deposition model the mass cancels,
%   so it only changes the reported marrow TIA (MBq h) and any cross-dose term
%   you compute from mass-scaled S-values.
arguments
    sex (1,1) string
end
switch upper(strtrim(sex))
    case "M",  m = 1170;
    case "F",  m = 900;
    otherwise
        error('rlt:badSex', 'Sex must be "M" or "F".');
end
end
