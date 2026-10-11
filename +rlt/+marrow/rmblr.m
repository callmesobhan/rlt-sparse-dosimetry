function [r, msg] = rmblr(hct, opts)
%RMBLR Red-marrow-to-blood (or -plasma) activity concentration ratio.
%   r = rlt.marrow.rmblr(0.47)                     % blood basis, RMECFF 0.19 -> 0.358
%   r = rlt.marrow.rmblr(0.45, RMECFF=0.19)
%   r = rlt.marrow.rmblr(Basis="plasma")           % plasma basis: r = RMECFF
%   [r, msg] = rlt.marrow.rmblr(hct)               % msg is "" unless Hct looks implausible
%
%   Blood basis:  RMBLR = RMECFF / (1 - Hct)   (Sgouros, J Nucl Med 1993)
%   Plasma basis: RMPR  = RMECFF
%   RMECFF is the red-marrow extracellular fluid fraction; 0.19 is the value
%   used in the literature for antibodies (with Hct 0.47 it gives RMBLR 0.36).
%   It was derived for activity that stays extracellular. Do not assume it holds
%   for small molecules with specific marrow uptake: see docs/MARROW.md.
arguments
    hct (1,1) double = NaN
    opts.RMECFF (1,1) double = 0.19
    opts.Basis (1,1) string {mustBeMember(opts.Basis, ["blood","plasma"])} = "blood"
end
if ~isfinite(opts.RMECFF) || opts.RMECFF <= 0 || opts.RMECFF >= 1
    error('rlt:badRMECFF', 'RMECFF must be strictly between 0 and 1.');
end
msg = "";
if opts.Basis == "plasma"
    r = opts.RMECFF;
    return
end
if ~isfinite(hct) || hct <= 0 || hct >= 1
    error('rlt:badHct', 'Hct must be a fraction strictly between 0 and 1 (0.45, not 45); got %g.', hct);
end
if hct < 0.25 || hct > 0.60
    msg = sprintf("Hct %.2f is outside the usual 0.25-0.60 range; please check.", hct);
end
r = opts.RMECFF / (1 - hct);
end
