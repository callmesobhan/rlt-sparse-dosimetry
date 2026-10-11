function s = lookupS(S, target, source, isotope)
%LOOKUPS Fetch one S-value (Gy per MBq h) from a table made by rlt.io.readSValues.
%   s = rlt.io.lookupS(S, "RedMarrow", "Kidneys", "Lu177")
%   Errors if the combination is missing, so a forgotten cross-dose term cannot
%   silently count as zero.
arguments
    S table
    target (1,1) string
    source (1,1) string
    isotope (1,1) string
end
isotope = rlt.isotope(isotope).name;
k = S.Target == target & S.Source == source & S.Isotope == isotope;
if ~any(k)
    error('rlt:missingSValue', 'No S-value for target "%s" <- source "%s" (%s) in the table.', ...
        target, source, isotope);
end
s = S.S_GyPerMBqh(find(k, 1));
end
