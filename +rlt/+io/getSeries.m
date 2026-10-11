function s = getSeries(T, patient, cycle, region)
%GETSERIES Extract one time-activity curve from a table made by rlt.io.readActivityTable.
%   s = rlt.io.getSeries(T, "P001", 1, "blood")
%
%   s.t        times (h), ascending
%   s.A        activity actually present at those times (decay correction undone)
%   s.unit     "MBq" | "MBq/mL" | "MBq/g"
%   s.isotope  isotope name
%   s.n        number of points
%   Passing s.t and s.A to rlt.fitTAC or rlt.marrow.selfDose is the intended path.
arguments
    T table
    patient (1,1) string
    cycle (1,1) double
    region (1,1) string
end
needed = ["PatientID","Cycle","Region","Time_h","ValuePhysical","UnitNorm","Isotope"];
if ~all(ismember(needed, string(T.Properties.VariableNames)))
    error('rlt:notNormalised', 'T must come from rlt.io.readActivityTable.');
end
k = T.PatientID == patient & T.Cycle == cycle & T.Region == region;
if ~any(k)
    avail = unique(T.Region(T.PatientID == patient & T.Cycle == cycle));
    if isempty(avail)
        error('rlt:seriesNotFound', 'No rows for patient "%s", cycle %d.', patient, cycle);
    end
    error('rlt:seriesNotFound', 'No region "%s" for patient "%s", cycle %d. Available: %s.', ...
        region, patient, cycle, strjoin(avail, ", "));
end
sub = sortrows(T(k, :), 'Time_h');
s.t       = sub.Time_h;
s.A       = sub.ValuePhysical;
s.unit    = sub.UnitNorm(1);
s.isotope = sub.Isotope(1);
s.n       = height(sub);
end
