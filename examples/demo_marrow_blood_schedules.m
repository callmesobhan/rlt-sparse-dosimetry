%% How many blood samples, and when, for blood-based marrow dosimetry?
% Scores blood sampling schedules on a virtual cohort whose true blood curve has
% three phases, so mono- and bi-exponential fits show model-form bias as well as
% noise. The cohort parameters are ILLUSTRATIVE (rlt.internal.presetParams):
% swap in distributions fitted to your own patients before trusting the ranking.
% Run from the repo root:  addpath(pwd); run("examples/demo_marrow_blood_schedules.m")
% Takes a few minutes because "select" fits up to three models per patient.

train = rlt.simulatePopulation(Preset="blood", N=100, Seed=99);
test  = rlt.simulatePopulation(Preset="blood", N=100, Seed=1);

schedules = { ...
    [1 24], ...
    [1 4 24], ...
    [0.5 2 6 24], ...
    [0.5 1 2 4 24 48], ...
    [0.5 1 2 4 8 24 48 72]};

% 5% per-sample CV is an assumption for a well-calibrated gamma counter.
R = rlt.evaluateSchedules(schedules, test, ...
    Methods=["mono","biexp","hybrid","select"], ...
    CV=0.05, LambdaPop=train.lambdaPop);

disp(R)
rlt.plotSchedules(R);
