function issues = validateActivityTable(T)
%VALIDATEACTIVITYTABLE Check a long-format activity table and list every problem at once.
%   issues = rlt.io.validateActivityTable(T) returns a table with columns
%   Severity ("error" | "warning"), DataRow, Message. An empty result means the table
%   is clean.
%
%   Required columns
%     PatientID       any text or number
%     Cycle           positive integer
%     Region          e.g. "blood", "wholebody", "vertebrae_L2L4", "kidney_L"
%     Time_h          hours after administration, >= 0
%     Value           >= 0
%     Unit            Bq, kBq, MBq, Bq/mL, kBq/mL, MBq/mL, Bq/g, kBq/g, MBq/g
%     DecayCorrected  1 if Value was decay-corrected to the administration time,
%                     0 if it is the activity actually present at Time_h. Required,
%                     never assumed: a wrong guess biases every TIA.
%   Optional column
%     Isotope         "Lu177", "Y90", "I131"
%
%   A (PatientID, Cycle, Region) group is one time-activity curve.
sev = strings(0, 1);
row = zeros(0, 1);
msg = strings(0, 1);

req   = ["PatientID","Cycle","Region","Time_h","Value","Unit","DecayCorrected"];
names = string(T.Properties.VariableNames);
missing = req(~ismember(req, names));
if ~isempty(missing)
    add("error", NaN, "Missing required column(s): " + strjoin(missing, ", ") + ".");
    issues = rlt.internal.issueTable(sev, row, msg);
    return
end
n = height(T);
if n == 0
    add("error", NaN, "The table has no rows.");
    issues = rlt.internal.issueTable(sev, row, msg);
    return
end

pid   = string(T.PatientID);
cyc   = double(T.Cycle);
reg   = string(T.Region);
time  = double(T.Time_h);
val   = double(T.Value);
unit  = string(T.Unit);
dc    = double(T.DecayCorrected);

for r = find(~(isfinite(cyc) & cyc >= 1 & cyc == round(cyc)))'
    add("error", r, "Cycle must be a positive integer.");
end
for r = find(ismissing(reg) | strlength(strtrim(reg)) == 0)'
    add("error", r, "Region is empty.");
end
for r = find(ismissing(pid) | strlength(strtrim(pid)) == 0)'
    add("error", r, "PatientID is empty.");
end
for r = find(~isfinite(time) | time < 0)'
    add("error", r, "Time_h must be a finite number >= 0 (hours after administration).");
end
for r = find(~isfinite(val) | val < 0)'
    add("error", r, "Value must be a finite number >= 0. Handle background-subtracted negatives explicitly before importing.");
end
for r = find(~ismember(dc, [0 1]))'
    add("error", r, "DecayCorrected must be 0 or 1.");
end

kind = strings(n, 1);
for r = 1:n
    info = rlt.internal.unitInfo(unit(r));
    if isempty(info)
        add("error", r, sprintf('Unknown unit "%s". Allowed: Bq, kBq, MBq, Bq/mL, kBq/mL, MBq/mL, Bq/g, kBq/g, MBq/g.', unit(r)));
    else
        kind(r) = info.kind;
    end
end

hasIso = ismember("Isotope", names);
if hasIso
    iso = string(T.Isotope);
    for r = 1:n
        try
            rlt.isotope(iso(r));
        catch
            add("error", r, sprintf('Unknown isotope "%s".', iso(r)));
        end
    end
end

key = pid + "|" + string(cyc) + "|" + reg + "|" + string(time);
[~, ~, ic] = unique(key);
cnt = accumarray(ic, 1);
for r = find(cnt(ic) > 1)'
    add("error", r, "Duplicate (PatientID, Cycle, Region, Time_h).");
end

gkey = pid + "|" + string(cyc) + "|" + reg;
for g = unique(gkey)'
    idx = find(gkey == g);
    label = strrep(g, "|", " / ");
    if numel(unique(kind(idx))) > 1
        add("error", idx(1), sprintf('Group "%s" mixes units of different kinds (activity vs per mL vs per g).', label));
    end
    if numel(unique(dc(idx))) > 1
        add("error", idx(1), sprintf('Group "%s" mixes decay-corrected and uncorrected rows.', label));
    end
    if hasIso && numel(unique(iso(idx))) > 1
        add("error", idx(1), sprintf('Group "%s" mixes isotopes.', label));
    end
    if numel(idx) < 2
        add("warning", idx(1), sprintf('Group "%s" has a single time point; most estimators need at least 2.', label));
    end
    if max(time(idx)) < 3
        add("warning", idx(1), sprintf('Group "%s": every time is below 3 h. Are times in hours (not days or minutes)?', label));
    end
    if ~issorted(time(idx))
        add("warning", idx(1), sprintf('Group "%s" is not in time order; it will be sorted.', label));
    end
end

issues = rlt.internal.issueTable(sev, row, msg);

    function add(s, r, m)
        sev(end+1, 1) = s;
        row(end+1, 1) = r;
        msg(end+1, 1) = m;
    end
end
