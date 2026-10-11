function [T, issues] = readActivityTable(file, opts)
%READACTIVITYTABLE Read, validate and normalise a long-format activity CSV.
%   [T, issues] = rlt.io.readActivityTable("activity.csv")
%   [T, issues] = rlt.io.readActivityTable("activity.csv", Isotope="Y90")
%
%   See rlt.io.validateActivityTable for the column definitions and
%   data/templates/activity_table_template.csv for an example. Any validation
%   error stops the import with the full list of problems; warnings are
%   reported through warning('rlt:tableWarning', ...) and returned in issues.
%
%   T is sorted by PatientID, Cycle, Region, Time_h and gains these columns:
%     Isotope        canonical name (from the file, or the Isotope option)
%     ValueNorm      Value converted to MBq, MBq/mL or MBq/g
%     UnitNorm       "MBq" | "MBq/mL" | "MBq/g"
%     ValuePhysical  ValueNorm with any decay correction undone, i.e. the
%                    activity actually present at Time_h. Fit this, not Value.
arguments
    file (1,1) string
    opts.Isotope (1,1) string = "Lu177"
end
if ~isfile(file)
    error('rlt:fileNotFound', 'File not found: %s', file);
end
raw    = rlt.internal.readCsvTable(file, ...
    ["PatientID","Region","Unit","Isotope"], ...
    ["Cycle","Time_h","Value","DecayCorrected"]);
issues = rlt.io.validateActivityTable(raw);

errs = issues(issues.Severity == "error", :);
if ~isempty(errs)
    error('rlt:invalidTable', 'Activity table "%s" has %d error(s):\n%s', ...
        file, height(errs), rlt.internal.formatIssues(errs));
end
warns = issues(issues.Severity == "warning", :);
if ~isempty(warns)
    warning('rlt:tableWarning', 'Activity table "%s" has %d warning(s):\n%s', ...
        file, height(warns), rlt.internal.formatIssues(warns));
end

n = height(raw);
T = raw;
T.PatientID      = string(raw.PatientID);
T.Region         = string(raw.Region);
T.Unit           = string(raw.Unit);
T.Cycle          = double(raw.Cycle);
T.Time_h         = double(raw.Time_h);
T.Value          = double(raw.Value);
T.DecayCorrected = logical(double(raw.DecayCorrected));
if ismember("Isotope", string(raw.Properties.VariableNames))
    iso = string(raw.Isotope);
    for r = 1:n
        iso(r) = rlt.isotope(iso(r)).name;
    end
else
    iso = repmat(rlt.isotope(opts.Isotope).name, n, 1);
end
T.Isotope = iso;

valueNorm = zeros(n, 1);
unitNorm  = strings(n, 1);
valuePhys = zeros(n, 1);
for r = 1:n
    info = rlt.internal.unitInfo(T.Unit(r));
    valueNorm(r) = T.Value(r) * info.factor;
    unitNorm(r)  = info.normUnit;
    if T.DecayCorrected(r)
        valuePhys(r) = rlt.undoDecayCorrection(T.Time_h(r), valueNorm(r), T.Isotope(r));
    else
        valuePhys(r) = valueNorm(r);
    end
end
T.ValueNorm     = valueNorm;
T.UnitNorm      = unitNorm;
T.ValuePhysical = valuePhys;

T = sortrows(T, {'PatientID', 'Cycle', 'Region', 'Time_h'});
end
