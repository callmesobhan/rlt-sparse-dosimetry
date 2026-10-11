function pop = simulatePopulation(opts)
%SIMULATEPOPULATION Virtual patients with known ground-truth time-integrated activity.
%   pop = rlt.simulatePopulation(Preset="kidney", N=200, Seed=1)
%
%   Presets: "kidney", "tumour" (washout plus an uptake term, so A(0) = 0) and
%   "blood" (three-phase washout from injection, for bone-marrow blood-based
%   dosimetry). Each patient has
%       A(t) = sum_k a_k*exp(-l_k*t)  [ - (sum_k a_k)*exp(-ku*t)  if uptake ]
%   with l_1 <= ... <= l_K and, for uptake presets, ku >= 1.2*l_K, so A(t) >= 0.
%   Parameters are drawn from the distributions in rlt.internal.presetParams,
%   which are ILLUSTRATIVE, not measured data.
%
%   Output struct
%     N, preset, isotope, units, nComp, hasUptake
%     amp          N x 1  total amplitude (sum of a_k)
%     a, lambda    N x K  true amplitudes and rates (1/h), slowest first
%     ku           N x 1  uptake rate (Inf when the preset has no uptake)
%     trueTIA      N x 1  exact integral over 0..Inf (units of A times h)
%     lambdaPop    population effective rate (1/h): median of amp/trueTIA, the
%                  rate of a mono-exponential with the same amplitude and TIA
%     eval         function handle: eval(tq) -> N x numel(tq) true activity
arguments
    opts.Preset (1,1) string = "kidney"
    opts.N (1,1) double {mustBeInteger, mustBePositive} = 200
    opts.Seed (1,1) double = 1
    opts.Isotope (1,1) string = "Lu177"
end
P   = rlt.internal.presetParams(opts.Preset);
iso = rlt.isotope(opts.Isotope);
N   = opts.N;
K   = P.nComp;
s   = RandStream('mt19937ar', 'Seed', opts.Seed);

amp = P.amp.median .* exp(log(P.amp.gsd) .* randn(s, N, 1));

w = P.frac .* exp(log(P.fracGsd) .* randn(s, N, K));
w = w ./ sum(w, 2);

Teff = P.Teff(:, 1)' .* exp(log(P.Teff(:, 2))' .* randn(s, N, K));
Teff = min(Teff, 0.98*iso.halfLifeH);              % never longer than physical
lam  = sort(log(2) ./ Teff, 2);                    % ascending: slowest first
a    = amp .* w;

if P.hasUptake
    Tup = P.Tup.median .* exp(log(P.Tup.gsd) .* randn(s, N, 1));
    ku  = max(log(2) ./ Tup, 1.2*lam(:, K));
    trueTIA = sum(a ./ lam, 2) - amp ./ ku;
else
    ku = inf(N, 1);
    trueTIA = sum(a ./ lam, 2);
end

pop.N         = N;
pop.preset    = opts.Preset;
pop.isotope   = iso.name;
pop.units     = P.units;
pop.nComp     = K;
pop.hasUptake = P.hasUptake;
pop.amp       = amp;
pop.a         = a;
pop.lambda    = lam;
pop.ku        = ku;
pop.trueTIA   = trueTIA;
pop.lambdaPop = median(amp ./ trueTIA);
if P.hasUptake
    pop.eval = @(tq) sumExp(a, lam, tq) - amp .* exp(-ku .* tq(:)');
else
    pop.eval = @(tq) sumExp(a, lam, tq);
end
end

function y = sumExp(a, lam, tq)
y = zeros(size(a, 1), numel(tq));
for k = 1:size(a, 2)
    y = y + a(:, k) .* exp(-lam(:, k) .* tq(:)');
end
end
