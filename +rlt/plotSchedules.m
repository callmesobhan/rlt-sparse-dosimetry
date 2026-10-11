function fig = plotSchedules(R)
%PLOTSCHEDULES Grouped bar charts of rlt.evaluateSchedules output.
%   fig = rlt.plotSchedules(R) shows, per schedule and estimator, the share of
%   virtual patients with |TIA error| <= 10% (top) and the 90th percentile of
%   |TIA error| (bottom).
arguments
    R table
end
schedules = unique(R.Schedule, 'stable');
methods   = unique(R.Method,   'stable');
within = nan(numel(schedules), numel(methods));
p90    = nan(numel(schedules), numel(methods));
for i = 1:numel(schedules)
    for j = 1:numel(methods)
        k = R.Schedule == schedules(i) & R.Method == methods(j);
        within(i, j) = R.Within10Pct(k);
        p90(i, j)    = R.P90AbsErrPct(k);
    end
end

fig = figure('Name', 'Sampling schedule comparison', 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');

ax = nexttile(tl);
bar(ax, within);
ylabel(ax, 'Patients within \pm10% (%)');
ylim(ax, [0 100]);
xticks(ax, 1:numel(schedules));
xticklabels(ax, cellstr(schedules));
legend(ax, cellstr(methods), 'Location', 'southeast');
grid(ax, 'on');

ax = nexttile(tl);
bar(ax, p90);
ylabel(ax, '90th pct |TIA error| (%)');
xticks(ax, 1:numel(schedules));
xticklabels(ax, cellstr(schedules));
xlabel(ax, 'Scan times after administration');
grid(ax, 'on');
end
