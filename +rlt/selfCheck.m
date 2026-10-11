function ok = selfCheck()
%SELFCHECK Quick installation and numerics check; needs no test framework.
%   ok = rlt.selfCheck() runs a handful of known-answer checks, prints PASS/FAIL
%   for each, and returns true only if all pass. Run it first on any new MATLAB
%   installation; the full suite is results = runtests("tests").
names  = strings(0, 1);
passed = false(0, 1);
detail = strings(0, 1);

% 1. MATLAB release
try
    old = isMATLABReleaseOlderThan("R2021a");
catch
    old = verLessThan('matlab', '9.10');
end
record("MATLAB R2021a or newer", ~old, string(version('-release')));

% 2. isotope constants
try
    iso = rlt.isotope("Lu177");
    record("Lu-177 decay constant = ln2 / half-life", ...
        abs(iso.lambdaPhys - log(2)/iso.halfLifeH) < 1e-15, "");
catch ME
    record("Lu-177 decay constant", false, string(ME.message));
end

% 3. dose constant
try
    d = rlt.absorbedDose(1, 1);
    record("1 MBq h in 1 kg of Lu-177 gives 7.729e-5 Gy", abs(d/7.7289e-5 - 1) < 1e-3, sprintf("%.6g Gy", d));
catch ME
    record("dose constant", false, string(ME.message));
end

% 4. mono-exponential recovery from two points
try
    lam = rlt.isotope("Lu177").lambdaPhys + 0.01;
    t = [24; 96];
    f = rlt.fitTAC(t, 1000*exp(-lam*t), Model="mono");
    record("mono fit recovers exact TIA from 2 points", abs(f.tia/(1000/lam) - 1) < 1e-3, sprintf("rel. error %.2g", f.tia/(1000/lam) - 1));
catch ME
    record("mono fit", false, string(ME.message));
end

% 5. bi-exponential recovery (noise-free)
try
    lp = rlt.isotope("Lu177").lambdaPhys;
    t = [1; 4; 8; 24; 48; 96; 168];
    A = 600*exp(-(lp+0.005)*t) + 400*exp(-(lp+0.1)*t);
    f = rlt.fitTAC(t, A, Model="biexp");
    truth = 600/(lp+0.005) + 400/(lp+0.1);
    record("biexp fit recovers noise-free TIA within 2%", abs(f.tia/truth - 1) < 0.02, sprintf("rel. error %.2g", f.tia/truth - 1));
catch ME
    record("biexp fit", false, string(ME.message));
end

% 6. RMBLR
try
    r = rlt.marrow.rmblr(0.47);
    record("RMBLR(Hct 0.47, RMECFF 0.19) = 0.358 (Sgouros: 0.36)", abs(r - 0.19/0.53) < 1e-12, sprintf("%.4f", r));
catch ME
    record("RMBLR", false, string(ME.message));
end

% 7. marrow self-dose does not depend on marrow mass
try
    t = [2; 24; 48; 96];
    C = 1e-3 * exp(-(rlt.isotope("Lu177").lambdaPhys + 0.02) * t);
    a = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", MarrowMassG=800);
    b = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", MarrowMassG=1500);
    record("marrow self-dose is independent of marrow mass", abs(a.doseGy/b.doseGy - 1) < 1e-12, "");
catch ME
    record("marrow mass independence", false, string(ME.message));
end

% 8. virtual-cohort truth against numerical integration
try
    pop = rlt.simulatePopulation(Preset="blood", N=3, Seed=2);
    tt = (0:0.25:6000)';
    num = trapz(tt, pop.eval(tt), 2);
    record("virtual-cohort TIA matches numerical integral", max(abs(num./pop.trueTIA - 1)) < 5e-3, "");
catch ME
    record("virtual cohort", false, string(ME.message));
end

% 9. data import round trip with decay-correction handling
try
    folder = tempname;
    mkdir(folder);
    cleaner = onCleanup(@() rmdir(folder, 's'));
    file = fullfile(folder, "a.csv");
    tbl = table(["P1";"P1"], [1;1], ["blood";"blood"], [4;24], [10;5], ["kBq/mL";"kBq/mL"], [1;1], ...
        'VariableNames', {'PatientID','Cycle','Region','Time_h','Value','Unit','DecayCorrected'});
    writetable(tbl, file);
    T = rlt.io.readActivityTable(file);
    s = rlt.io.getSeries(T, "P1", 1, "blood");
    expected = 5e-3 * exp(-rlt.isotope("Lu177").lambdaPhys * 24);
    record("import converts units and undoes decay correction", abs(s.A(2)/expected - 1) < 1e-12 && s.unit == "MBq/mL", "");
catch ME
    record("data import", false, string(ME.message));
end

ok = all(passed);
fprintf("\nrlt.selfCheck: %d/%d passed (rlt %s)\n", sum(passed), numel(passed), rlt.version());

    function record(name, cond, info)
        names(end+1, 1)  = name;
        passed(end+1, 1) = logical(cond);
        detail(end+1, 1) = info;
        if cond
            tag = "PASS";
        else
            tag = "FAIL";
        end
        if strlength(info) > 0
            fprintf("  [%s] %s  (%s)\n", tag, name, info);
        else
            fprintf("  [%s] %s\n", tag, name);
        end
    end
end
