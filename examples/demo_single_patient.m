%% Single-patient workflow: sparse scans -> model -> TIA with uncertainty -> dose
% Uses a virtual patient so the ground truth is known. Replace t and A with
% your own quantified scan data (activity in MBq) to analyse a real patient.
% Run from the repo root:  addpath(pwd); run("examples/demo_single_patient.m")

isotope  = "Lu177";
injected = 7400;                       % MBq administered
massKg   = 0.30;                       % reference kidney mass; use patient-specific volume

pop   = rlt.simulatePopulation(Preset="kidney", N=1, Seed=11);
t     = [4 24 72 168];                 % scan times after administration (h)
truth = injected * pop.eval(t);        % MBq in the organ
rng(1);
A     = truth(:) .* exp(0.08*randn(numel(t),1));   % add ~8% quantification noise

%% Candidate models, and which one the data can justify
[best, cand] = rlt.selectModel(t, A, Isotope=isotope);
disp(cand)
fprintf("Selected model: %s\n", best.model);

%% TIA by several estimators, against the known truth
trueTIA = injected * pop.trueTIA;
methods = ["mono","hybrid","select"];
fprintf("\n%-8s %12s %10s\n", "Method", "TIA (MBq h)", "error (%)");
for m = methods
    v = rlt.estimate(m, t, A, Isotope=isotope);
    fprintf("%-8s %12.0f %10.1f\n", m, v, 100*(v/trueTIA - 1));
end
fprintf("%-8s %12.0f\n", "truth", trueTIA);

%% Uncertainty from 10% per-scan quantification error
u = rlt.propagateTIA(t, A, Method="select", CV=0.10, N=500, Isotope=isotope);
fprintf("\nTIA median %.0f MBq h, 95%% interval [%.0f, %.0f] (+/- %.0f%%)\n", ...
    u.median, u.ci95(1), u.ci95(2), 100*u.relHalfWidth);

%% Absorbed dose
D = rlt.absorbedDose(best.tia, massKg, Isotope=isotope);
Dlo = rlt.absorbedDose(u.ci95(1), massKg, Isotope=isotope);
Dhi = rlt.absorbedDose(u.ci95(2), massKg, Isotope=isotope);
fprintf("Absorbed dose %.2f Gy (95%% interval %.2f to %.2f Gy), local-deposition model\n", D, Dlo, Dhi);

%% Plot
tq = linspace(0, 400, 400);
figure('Color', 'w'); hold on
plot(tq, injected*pop.eval(tq), 'k-', 'DisplayName', 'truth');
plot(tq, best.predict(tq), 'b--', 'DisplayName', "fit: " + best.model);
plot(t, A, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'measured');
xlabel('Time after administration (h)'); ylabel('Activity in organ (MBq)');
legend('Location', 'northeast'); grid on
