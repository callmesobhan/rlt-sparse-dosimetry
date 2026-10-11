function [S, issues] = readSValues(file)
%READSVALUES Read and validate a table of S-values you obtained from a trusted source.
%   [S, issues] = rlt.io.readSValues("svalues.csv")
%
%   The toolbox ships NO S-values: they depend on phantom, isotope and source
%   model, and must come from OLINDA/EXM, MIRDcalc, IDAC-Dose or your own Monte
%   Carlo, with the reference recorded.
%
%   Required columns
%     Target      region receiving the dose, e.g. "RedMarrow"
%     Source      region holding the activity, e.g. "RedMarrow", "Kidneys", "RemainderOfBody"
%     Isotope     "Lu177", "Y90", "I131"
%     S           the value
%     Unit        "Gy/(MBq h)", "mGy/(MBq h)", "Gy/(MBq s)" or "mGy/(MBq s)"
%     Reference   where the value came from (software, version, phantom); must not be blank
%   Optional column
%     ReferenceMass_g   mass of the Source (or Target) in the phantom, for mass scaling
%
%   S gains the column S_GyPerMBqh. Errors are listed together and stop the import.
arguments
    file (1,1) string
end
if ~isfile(file)
    error('rlt:fileNotFound', 'File not found: %s', file);
end
raw   = rlt.internal.readCsvTable(file, ...
    ["Target","Source","Isotope","Unit","Reference"], ...
    ["S","ReferenceMass_g"]);
names = string(raw.Properties.VariableNames);

req = ["Target","Source","Isotope","S","Unit","Reference"];
missing = req(~ismember(req, names));
if ~isempty(missing)
    error('rlt:invalidTable', 'S-value table "%s" is missing column(s): %s.', file, strjoin(missing, ", "));
end

unitNames = ["Gy/(MBq h)", "mGy/(MBq h)", "Gy/(MBq s)", "mGy/(MBq s)"];
unitScale = [1, 1e-3, 3600, 3.6];       % to Gy/(MBq h)

n = height(raw);
sev = strings(0, 1);  row = zeros(0, 1);  msg = strings(0, 1);
if n == 0
    sev = "error"; row = NaN; msg = "The table has no rows.";
end
val  = double(raw.S);
unit = string(raw.Unit);
ref  = string(raw.Reference);
iso  = string(raw.Isotope);
sGy  = nan(n, 1);
for r = 1:n
    k = find(unitNames == unit(r), 1);
    if isempty(k)
        sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
        msg(end+1, 1) = sprintf('Unknown unit "%s". Allowed: %s.', unit(r), strjoin(unitNames, ", "));
    elseif ~isfinite(val(r)) || val(r) < 0
        sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
        msg(end+1, 1) = "S must be a finite number >= 0.";
    else
        sGy(r) = val(r) * unitScale(k);
    end
    if ismissing(ref(r)) || strlength(strtrim(ref(r))) == 0
        sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
        msg(end+1, 1) = "Reference is blank: record where this S-value came from.";
    end
    try
        iso(r) = rlt.isotope(iso(r)).name;
    catch
        sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
        msg(end+1, 1) = sprintf('Unknown isotope "%s".', iso(r));
    end
end
key = string(raw.Target) + "|" + string(raw.Source) + "|" + iso;
[~, ~, ic] = unique(key);
cnt = accumarray(ic, 1);
for r = find(cnt(ic) > 1)'
    sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
    msg(end+1, 1) = "Duplicate (Target, Source, Isotope).";
end
if ismember("ReferenceMass_g", names)
    m = double(raw.ReferenceMass_g);
    for r = find(isinf(m) | m <= 0)'                 % blank (NaN) is allowed
        sev(end+1, 1) = "error"; row(end+1, 1) = r; %#ok<AGROW>
        msg(end+1, 1) = "ReferenceMass_g must be a positive number (or blank).";
    end
end

issues = rlt.internal.issueTable(sev, row, msg);
errs = issues(issues.Severity == "error", :);
if ~isempty(errs)
    error('rlt:invalidTable', 'S-value table "%s" has %d error(s):\n%s', ...
        file, height(errs), rlt.internal.formatIssues(errs));
end

isPlaceholder = contains(lower(ref), "placeholder");
if any(isPlaceholder)
    warning('rlt:placeholderSValue', ...
        '%d S-value row(s) are marked PLACEHOLDER in Reference (rows: %s). Replace them with values from a real source before using any result.', ...
        sum(isPlaceholder), strjoin(string(find(isPlaceholder)'), ", "));
end

S = raw;
S.Target  = string(raw.Target);
S.Source  = string(raw.Source);
S.Isotope = iso;
S.S_GyPerMBqh = sGy;
end
