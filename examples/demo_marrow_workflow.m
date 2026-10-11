%% Red-marrow dosimetry workflow, end to end, on the shipped SYNTHETIC template data
% The template patient is made up: no number below is a clinical result.
% Run from the repo root:  addpath(pwd); run("examples/demo_marrow_workflow.m")

root = pwd;
tdir = fullfile(root, "data", "templates");

%% 1. Import with validation (units, decay-correction flag, duplicates, ...)
T = rlt.io.readActivityTable(fullfile(tdir, "activity_table_template.csv"));
P = rlt.io.readPatientInfo(fullfile(tdir, "patient_info_template.csv"));
pid = "EXAMPLE_SYNTHETIC";
pat = P(P.PatientID == pid, :);
injected = 7400;                                    % MBq administered (synthetic)

blood = rlt.io.getSeries(T, pid, 1, "blood");       % decay correction already handled
fprintf("Blood: %d samples from %g to %g h, unit %s\n", blood.n, blood.t(1), blood.t(end), blood.unit);

%% 2. Blood-based self-dose, with the model comparison that matters
res = rlt.marrow.selfDose(blood.t, blood.A, ...
    ConcentrationUnit=blood.unit, Hct=pat.Hct, Sex=pat.Sex, ...
    AdministeredMBq=injected, Isotope=blood.isotope);
fprintf("\nSelected model: %s\n", res.fit.model);
fprintf("Self-dose to red marrow: %.3f Gy  (%.1f mGy/GBq)\n", res.doseGy, 1000*res.dosePerGBq);
disp(res.candidates)
fprintf("Model-to-model dose spread: %.0f%%\n", 100*res.modelSpread);

%% 3. Warnings and assumptions travel with the number
fprintf("\nWarnings:\n");  disp(res.warnings)
fprintf("Assumptions:\n"); disp(res.assumptions)

%% 4. What if the marrow-to-blood ratio rises with time? (what-if, not data)
rising = rlt.marrow.selfDose(blood.t, blood.A, ConcentrationUnit=blood.unit, Isotope=blood.isotope, ...
    RatioFcn=@(t) res.ratio .* (1 + 0.5*(1 - exp(-t/72))), AdministeredMBq=injected);
fprintf("Constant ratio: %.1f mGy/GBq; ratio rising by up to 50%%: %.1f mGy/GBq\n", ...
    1000*res.dosePerGBq, 1000*rising.dosePerGBq);

%% 5. Which input does the dose depend on?
fprintf("\nElasticity d ln D / d ln p:  Hct %.2f,  RMECFF %.2f\n", ...
    rlt.marrow.elasticity(blood.t, blood.A, "Hct", pat.Hct, "ConcentrationUnit", blood.unit), ...
    rlt.marrow.elasticity(blood.t, blood.A, "RMECFF", 0.19, "Hct", pat.Hct, "ConcentrationUnit", blood.unit));

%% 6. Add the other contributions. S-VALUES BELOW ARE TOY NUMBERS, NOT REAL.
% Use rlt.io.readSValues on values from OLINDA/EXM, MIRDcalc or IDAC-Dose.
wb = rlt.io.getSeries(T, pid, 1, "wholebody");
wbFit = rlt.fitTAC(wb.t, wb.A, Model="biexp", Isotope=wb.isotope);
rob = rlt.marrow.remainderOfBody(wbFit.tia, 0);     % no separately modelled organs in this demo
contrib = table(["Marrow self (blood-based)"; "Remainder of body"], ...
    [res.marrowTIA; rob], ...
    [res.doseGy/res.marrowTIA; 1e-7], ...           % second value: TOY S-value
    [0.25; 0.25], [0.20; 0.30], ...
    'VariableNames', {'Source','TIA_MBqh','S_GyPerMBqh','TIA_CV','S_CV'});
tot = rlt.marrow.totalDose(contrib);
unc = rlt.marrow.uncertainty(contrib, SharedTIACV=0.10);
disp(tot.table)
fprintf("Total: %.3f Gy, 95%% interval [%.3f, %.3f] Gy (toy cross-dose S-value)\n", ...
    tot.totalGy, unc.ci95(1), unc.ci95(2));
disp(unc.perSource)

%% 7. Auditable report
out = fullfile(tempdir, "rlt_marrow_report_demo.md");
rlt.report.marrow(out, res, Total=tot, Uncertainty=unc, Patient=pid);
fprintf("\nReport written to %s\n", out);
