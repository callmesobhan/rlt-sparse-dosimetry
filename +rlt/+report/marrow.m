function marrow(file, selfRes, opts)
%MARROW Write a plain-text (Markdown) marrow dosimetry report with full provenance.
%   rlt.report.marrow("report.md", selfRes, Total=totRes, Uncertainty=unc, Patient="P001")
%
%   selfRes      result of rlt.marrow.selfDose
%   Total        result of rlt.marrow.totalDose (optional)
%   Uncertainty  result of rlt.marrow.uncertainty (optional)
%   Patient      free-text label; do not put identifiable information in it
%                unless your data-protection rules allow it
%
%   The report lists inputs, every assumption and warning, the model comparison,
%   results and software provenance, so the number can be audited later. It is
%   not a clinical document.
arguments
    file (1,1) string
    selfRes (1,1) struct
    opts.Total struct = struct()
    opts.Uncertainty struct = struct()
    opts.Patient (1,1) string = "(not given)"
end
s = selfRes;
p = s.provenance;

L = strings(0, 1);
L(end+1) = "# Red-marrow dosimetry report";
L(end+1) = "";
L(end+1) = "Research use only. Not a clinical document.";
L(end+1) = "";
L(end+1) = "- Patient label: " + opts.Patient;
L(end+1) = "- Isotope: " + s.isotope;
L(end+1) = "- Method: " + s.source + "-based self-dose";
L(end+1) = "- Toolbox: " + p.toolbox + " " + p.version + ", MATLAB " + p.release + " (" + p.computer + "), " + p.timestamp;
L(end+1) = "";
L(end+1) = "## Result";
L(end+1) = "";
L(end+1) = sprintf("- Self-dose to red marrow: **%.4g Gy** (model: %s)", s.doseGy, s.fit.model);
if ~isnan(s.dosePerGBq)
    L(end+1) = sprintf("- Per administered activity: %.4g Gy/GBq", s.dosePerGBq);
end
if ~isnan(s.marrowTIA)
    L(end+1) = sprintf("- Marrow TIA: %.4g MBq h (marrow mass %.0f g)", s.marrowTIA, s.massG);
end
L(end+1) = sprintf("- Share of TIA extrapolated beyond the sampled window: %.0f%%", 100*s.fit.extrapolatedFraction);
if ~isnan(s.modelSpread)
    L(end+1) = sprintf("- Model-to-model dose spread: %.0f%%", 100*s.modelSpread);
end
L(end+1) = "";

if isfield(opts.Total, 'totalGy')
    L(end+1) = "## Total marrow dose (all listed sources)";
    L(end+1) = "";
    L(end+1) = sprintf("**Total: %.4g Gy**", opts.Total.totalGy);
    L(end+1) = "";
    L(end+1) = "| Source | TIA (MBq h) | S (Gy/(MBq h)) | Dose (Gy) | Share (%) |";
    L(end+1) = "|---|---|---|---|---|";
    T = opts.Total.table;
    for i = 1:height(T)
        L(end+1) = sprintf("| %s | %.4g | %.4g | %.4g | %.1f |", ...
            T.Source(i), T.TIA_MBqh(i), T.S_GyPerMBqh(i), T.DoseGy(i), T.SharePct(i));
    end
    L(end+1) = "";
end
if isfield(opts.Uncertainty, 'ci95')
    u = opts.Uncertainty;
    L(end+1) = "## Uncertainty (Monte Carlo)";
    L(end+1) = "";
    L(end+1) = sprintf("Median %.4g Gy, 95%% interval [%.4g, %.4g] Gy (half-width %.0f%% of median).", ...
        u.median, u.ci95(1), u.ci95(2), 100*u.relHalfWidth);
    L(end+1) = "";
    L(end+1) = "| Source | Median (Gy) | 95% interval (Gy) | Variance share (%) |";
    L(end+1) = "|---|---|---|---|";
    for i = 1:height(u.perSource)
        L(end+1) = sprintf("| %s | %.4g | [%.4g, %.4g] | %.1f |", u.perSource.Source(i), ...
            u.perSource.MedianGy(i), u.perSource.Lo95Gy(i), u.perSource.Hi95Gy(i), ...
            u.perSource.VarianceSharePct(i));
    end
    L(end+1) = "";
end

L(end+1) = "## Model comparison";
L(end+1) = "";
L(end+1) = "| Model | Parameters | Dose (Gy) | AICc | Selected |";
L(end+1) = "|---|---|---|---|---|";
c = s.candidates;
for i = 1:height(c)
    L(end+1) = sprintf("| %s | %d | %.4g | %.4g | %s |", c.Model(i), c.nParams(i), ...
        c.DoseGy(i), c.AICc(i), string(c.Selected(i)));
end
L(end+1) = "";

L(end+1) = "## Warnings";
L(end+1) = "";
if isempty(s.warnings)
    L(end+1) = "None.";
else
    for w = s.warnings'
        L(end+1) = "- " + w;
    end
end
L(end+1) = "";
L(end+1) = "## Assumptions";
L(end+1) = "";
for a = s.assumptions'
    L(end+1) = "- " + a;
end
L(end+1) = "";

fid = fopen(file, 'w');
if fid < 0
    error('rlt:cannotWrite', 'Cannot open "%s" for writing.', file);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', L);
end
