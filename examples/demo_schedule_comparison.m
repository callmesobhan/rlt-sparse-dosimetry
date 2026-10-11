%% Which imaging schedule and estimator keeps TIA error under control?
% Compares candidate post-therapy schedules on a virtual cohort with known truth.
% The cohort distributions are ILLUSTRATIVE presets: edit rlt.internal.presetParams
% (or build your own pop struct) with values from your own patients before you
% trust any conclusion for a real protocol.
% Run from the repo root:  addpath(pwd); run("examples/demo_schedule_comparison.m")

preset = "kidney";                      % or "tumour"

% Population rate for the single-time-point method comes from a DIFFERENT
% cohort than the one being scored, as it would in practice.
train = rlt.simulatePopulation(Preset=preset, N=300, Seed=99);
test  = rlt.simulatePopulation(Preset=preset, N=300, Seed=1);

schedules = { ...
    96, ...
    [24 96], ...
    [24 168], ...
    [4 24 96], ...
    [24 72 168], ...
    [4 24 48 96 168]};

R = rlt.evaluateSchedules(schedules, test, ...
    Methods=["mono","hybrid","stp","select"], ...
    CV=0.10, LambdaPop=train.lambdaPop);

disp(R)
fprintf("Population effective rate used for 'stp': %.4f 1/h -> optimal single scan at %.0f h\n", ...
    train.lambdaPop, rlt.optimalSinglePointTime(train.lambdaPop));

rlt.plotSchedules(R);
