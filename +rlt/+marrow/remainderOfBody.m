function [tiaROB, info] = remainderOfBody(tiaWB, tiaSources, opts)
%REMAINDEROFBODY Time-integrated activity of the remainder of the body.
%   tiaROB = rlt.marrow.remainderOfBody(tiaWB, [tiaKidneys tiaLiver tiaTumours ...])
%
%   Remainder of body = whole-body TIA minus the TIA of every region treated as
%   a separate source (the marrow-dosimetry convention: each source region is
%   counted once). Fit the whole-body curve with a decreasing model (rlt.fitTAC,
%   Model="mono" or "biexp") on NOT decay-corrected whole-body activity, and
%   make sure every term uses the same units (MBq h).
%
%   Negative Option: if the sources add up to more than the whole body, the inputs
%   are inconsistent (calibration, units, decay correction, or an unmodelled
%   source). Default "error" stops; "clip" sets the remainder to 0 and reports it.
%
%   info fields: sumSources, fractionOfWholeBody, clipped.
arguments
    tiaWB (1,1) double {mustBeNonnegative, mustBeFinite}
    tiaSources double {mustBeNonnegative, mustBeFinite}
    opts.Negative (1,1) string {mustBeMember(opts.Negative, ["error","clip"])} = "error"
end
tot = sum(tiaSources(:));
tiaROB = tiaWB - tot;
info.sumSources = tot;
if tiaWB > 0
    info.fractionOfWholeBody = tot / tiaWB;
else
    info.fractionOfWholeBody = NaN;
end
info.clipped = false;
if tiaROB < 0
    if tiaROB > -1e-9 * max(tiaWB, 1)
        tiaROB = 0;
    elseif opts.Negative == "error"
        error('rlt:negativeROB', ...
            ['Source regions (%.4g MBq h) exceed the whole body (%.4g MBq h). ', ...
             'Check units, calibration, decay correction and the whole-body model.'], tot, tiaWB);
    else
        tiaROB = 0;
        info.clipped = true;
    end
end
end
