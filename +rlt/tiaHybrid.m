function tia = tiaHybrid(t, A, opts)
%TIAHYBRID Time-integrated activity by trapezoidal rule plus an analytic tail.
%   tia = rlt.tiaHybrid(t, A)
%
%   The measured part is integrated with the trapezoidal rule. The curve before
%   the first scan is closed according to Origin, and the part after the last scan
%   with an exponential tail.
%
%   Name-value options
%     Isotope : "Lu177" (default), "Y90", "I131"
%     Origin  : how to close the curve before the first sample, if t(1) > 0
%               "zero"  (default) adds (0, 0), i.e. activity builds up linearly.
%                       Suitable for organs and lesions that fill after injection.
%               "flat"  adds (0, A(1)), i.e. instantaneous uptake.
%               "decay" extends the first two samples back to t = 0 along the
%                       exponential through them (falls back to "flat" if the
%                       curve is not falling). Use it for blood and whole body,
%                       whose activity is highest at injection; "zero" would
%                       badly underestimate those.
%     Tail    : "effective" (default) uses the washout rate between the last two
%               scans, floored at the physical rate; "physical" uses physical
%               decay only (an upper bound on the clearance rate).
arguments
    t double {mustBeVector, mustBeNonnegative, mustBeFinite}
    A double {mustBeVector, mustBeNonnegative, mustBeFinite}
    opts.Isotope (1,1) string = "Lu177"
    opts.Origin (1,1) string {mustBeMember(opts.Origin, ["zero","flat","decay"])} = "zero"
    opts.Tail (1,1) string {mustBeMember(opts.Tail, ["effective","physical"])} = "effective"
end
t = t(:);
A = A(:);
rlt.internal.checkSeries(t, A);
iso   = rlt.isotope(opts.Isotope);
nMeas = numel(t);

head = 0;                               % integral over [0, t(1)]
if t(1) > 0
    origin = opts.Origin;
    if origin == "decay"
        if A(2) < A(1) && A(2) > 0
            lam0 = log(A(1)/A(2)) / (t(2) - t(1));
            head = A(1) * (exp(lam0*t(1)) - 1) / lam0;
        else
            origin = "flat";
        end
    end
    switch origin
        case "zero"
            head = 0.5 * t(1) * A(1);
        case "flat"
            head = t(1) * A(1);
    end
end

body    = head + trapz(t, A);
lamTail = iso.lambdaPhys;
if opts.Tail == "effective" && nMeas >= 2
    if A(nMeas) > 0 && A(nMeas-1) > A(nMeas)
        lamEff  = log(A(nMeas-1)/A(nMeas)) / (t(nMeas) - t(nMeas-1));
        lamTail = max(lamEff, iso.lambdaPhys);
    end
end
tia = body + A(nMeas)/lamTail;
end
