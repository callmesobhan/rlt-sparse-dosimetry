function [P, issues] = readPatientInfo(file)
%READPATIENTINFO Read and validate per-patient covariates (haematocrit, sex, mass, ...).
%   [P, issues] = rlt.io.readPatientInfo("patients.csv")
%
%   Required columns
%     PatientID   unique per row
%     Sex         "M" or "F" (selects the reference red-marrow mass)
%     Hct         haematocrit as a FRACTION (0.45, not 45)
%     Mass_kg     body mass in kg
%   Optional columns
%     Height_cm
%     BoneMetastases   yes | no | unknown (or 1 | 0). Normalised to text. Needed
%                      to decide whether blood-based or image-based marrow
%                      dosimetry is trustworthy.
%
%   Errors stop the import and are listed together; warnings (implausible but
%   possible values) are returned in issues and raised as warning('rlt:patientWarning').
arguments
    file (1,1) string
end
if ~isfile(file)
    error('rlt:fileNotFound', 'File not found: %s', file);
end
raw   = readtable(file, 'TextType', 'string');
names = string(raw.Properties.VariableNames);

sev = strings(0, 1);
row = zeros(0, 1);
msg = strings(0, 1);

req = ["PatientID","Sex","Hct","Mass_kg"];
missing = req(~ismember(req, names));
if ~isempty(missing)
    add("error", NaN, "Missing required column(s): " + strjoin(missing, ", ") + ".");
    finish();
end

n    = height(raw);
pid  = string(raw.PatientID);
sex  = upper(strtrim(string(raw.Sex)));
hct  = double(raw.Hct);
mass = double(raw.Mass_kg);

if n == 0
    add("error", NaN, "The table has no rows.");
    finish();
end
[~, ~, ic] = unique(pid);
cnt = accumarray(ic, 1);
for r = find(cnt(ic) > 1)'
    add("error", r, "Duplicate PatientID.");
end
for r = find(~ismember(sex, ["M","F"]))'
    add("error", r, 'Sex must be "M" or "F".');
end
for r = find(~isfinite(hct) | hct <= 0 | hct >= 1)'
    add("error", r, "Hct must be a fraction strictly between 0 and 1 (0.45, not 45).");
end
for r = find(isfinite(hct) & hct > 0 & hct < 1 & (hct < 0.25 | hct > 0.60))'
    add("warning", r, sprintf("Hct %.2f is outside the usual 0.25-0.60 range; please check.", hct(r)));
end
for r = find(~isfinite(mass) | mass <= 0)'
    add("error", r, "Mass_kg must be a positive number.");
end
for r = find(isfinite(mass) & mass > 0 & (mass < 30 | mass > 200))'
    add("warning", r, sprintf("Mass %.1f kg is outside 30-200 kg; please check.", mass(r)));
end
if ismember("Height_cm", names)
    h = double(raw.Height_cm);
    for r = find(~isfinite(h) | h < 100 | h > 230)'
        add("warning", r, "Height_cm is missing or outside 100-230 cm; please check.");
    end
end

bm = repmat("unknown", n, 1);
if ismember("BoneMetastases", names)
    b = lower(strtrim(string(raw.BoneMetastases)));
    for r = 1:n
        switch b(r)
            case {"1","yes","true"},      bm(r) = "yes";
            case {"0","no","false"},      bm(r) = "no";
            case {"unknown",""},          bm(r) = "unknown";
            otherwise
                add("error", r, 'BoneMetastases must be yes, no or unknown (or 1/0).');
        end
    end
end

finish();
P = raw;
P.PatientID = pid;
P.Sex       = sex;
P.Hct       = hct;
P.Mass_kg   = mass;
P.BoneMetastases = bm;
issues = rlt.internal.issueTable(sev, row, msg);

    function add(s, r, m)
        sev(end+1, 1) = s;
        row(end+1, 1) = r;
        msg(end+1, 1) = m;
    end

    function finish()
        iss  = rlt.internal.issueTable(sev, row, msg);
        errs = iss(iss.Severity == "error", :);
        if ~isempty(errs)
            error('rlt:invalidTable', 'Patient table "%s" has %d error(s):\n%s', ...
                file, height(errs), rlt.internal.formatIssues(errs));
        end
        warns = iss(iss.Severity == "warning", :);
        if ~isempty(warns)
            warning('rlt:patientWarning', 'Patient table "%s" has %d warning(s):\n%s', ...
                file, height(warns), rlt.internal.formatIssues(warns));
        end
    end
end
