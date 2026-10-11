function iso = isotope(name)
%ISOTOPE Physical constants for therapeutic radionuclides.
%   iso = rlt.isotope("Lu177") returns a struct with fields
%     name         canonical isotope name (string)
%     halfLifeH    physical half-life (h)
%     lambdaPhys   physical decay constant (1/h)
%     meanBetaMeV  mean beta-minus energy emitted per decay (MeV)
%
%   Supported: "Lu177", "Y90", "I131".
%
%   These are reference values for development. Verify against ICRP
%   Publication 107 / MIRD data before using results clinically.
arguments
    name (1,1) string
end

switch upper(name)
    case "LU177"
        canonical = "Lu177";  halfLifeH = 6.647*24;   meanBeta = 0.134;
    case "Y90"
        canonical = "Y90";    halfLifeH = 2.671*24;   meanBeta = 0.9337;
    case "I131"
        canonical = "I131";   halfLifeH = 8.021*24;   meanBeta = 0.182;
    otherwise
        error('rlt:unknownIsotope', ...
            'Unknown isotope "%s". Supported: Lu177, Y90, I131.', name);
end

iso = struct('name', canonical, ...
             'halfLifeH', halfLifeH, ...
             'lambdaPhys', log(2)/halfLifeH, ...
             'meanBetaMeV', meanBeta);
end
