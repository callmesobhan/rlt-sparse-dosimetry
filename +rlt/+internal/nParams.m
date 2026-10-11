function p = nParams(model)
%NPARAMS Number of free parameters (amplitudes + rates) of a kinetic model.
%   mono   : A(t) = a*exp(-l1*t)                              -> 2
%   uptake : A(t) = a*(exp(-l1*t) - exp(-l2*t))               -> 3
%   biexp  : A(t) = a1*exp(-l1*t) + a2*exp(-l2*t)             -> 4
%   triexp : A(t) = a1*exp(-l1*t) + a2*exp(-l2*t) + a3*exp(-l3*t) -> 6
switch string(model)
    case "mono",   p = 2;
    case "uptake", p = 3;
    case "biexp",  p = 4;
    case "triexp", p = 6;
    otherwise
        error('rlt:unknownModel', 'Unknown model "%s".', string(model));
end
end
