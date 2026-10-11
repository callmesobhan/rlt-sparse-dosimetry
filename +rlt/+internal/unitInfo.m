function info = unitInfo(unit)
%UNITINFO Conversion of an activity unit string to MBq-based units.
%   info = rlt.internal.unitInfo("kBq/mL") returns a struct with
%     unit      the matched unit string
%     factor    multiply a value in this unit by factor to get the normalised unit
%     kind      "activity" | "perMilliliter" | "perGram"
%     normUnit  "MBq" | "MBq/mL" | "MBq/g"
%   An unknown unit returns [].
names  = ["Bq","kBq","MBq","Bq/mL","kBq/mL","MBq/mL","Bq/g","kBq/g","MBq/g"];
factor = [1e-6, 1e-3, 1, 1e-6, 1e-3, 1, 1e-6, 1e-3, 1];
kind   = ["activity","activity","activity", ...
          "perMilliliter","perMilliliter","perMilliliter", ...
          "perGram","perGram","perGram"];
norm   = ["MBq","MBq","MBq","MBq/mL","MBq/mL","MBq/mL","MBq/g","MBq/g","MBq/g"];

k = find(names == string(unit), 1);
if isempty(k)
    info = [];
else
    info = struct('unit', names(k), 'factor', factor(k), ...
                  'kind', kind(k), 'normUnit', norm(k));
end
end
