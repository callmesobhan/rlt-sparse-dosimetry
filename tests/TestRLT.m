classdef TestRLT < matlab.unittest.TestCase
    %TESTRLT Core fitting, estimators, dose and virtual-cohort tests.
    %   Run all tests from the repo root with:  results = runtests("tests")

    methods (TestClassSetup)
        function addRepoToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)

        % ---- constants -----------------------------------------------------

        function isotopeConstants(tc)
            iso = rlt.isotope("lu177");
            tc.verifyEqual(iso.name, "Lu177");
            tc.verifyEqual(iso.halfLifeH, 6.647*24, 'RelTol', 1e-12);
            tc.verifyEqual(iso.lambdaPhys, log(2)/iso.halfLifeH, 'RelTol', 1e-12);
            tc.verifyError(@() rlt.isotope("Xx1"), 'rlt:unknownIsotope');
        end

        function versionIsAString(tc)
            tc.verifyClass(rlt.version(), 'string');
            tc.verifyTrue(rlt.selfCheck());
        end

        % ---- fitting -------------------------------------------------------

        function monoFitRecoversExactTIAFromTwoPoints(tc)
            iso = rlt.isotope("Lu177");
            lam = iso.lambdaPhys + 0.01;
            t = [24; 96];
            A = 1000 * exp(-lam * t);
            f = rlt.fitTAC(t, A, Model="mono");
            tc.verifyEqual(f.tia, 1000/lam, 'RelTol', 1e-3);
            tc.verifyTrue(f.exactlyDetermined);
            tc.verifyGreaterThanOrEqual(f.lambda, iso.lambdaPhys);
        end

        function rowAndColumnInputsGiveTheSameFit(tc)
            t = [24 48 96];
            A = [100 70 40];
            tc.verifyEqual(rlt.fitTAC(t, A).tia, rlt.fitTAC(t(:), A(:)).tia, 'RelTol', 1e-12);
        end

        function biexpFitRecoversNoiseFreeTIA(tc)
            iso = rlt.isotope("Lu177");
            l1 = iso.lambdaPhys + 0.005;
            l2 = iso.lambdaPhys + 0.100;
            t = [1; 4; 8; 24; 48; 96; 168];
            A = 600*exp(-l1*t) + 400*exp(-l2*t);
            f = rlt.fitTAC(t, A, Model="biexp");
            tc.verifyEqual(f.tia, 600/l1 + 400/l2, 'RelTol', 2e-2);
        end

        function triexpFitRecoversNoiseFreeTIA(tc)
            lp = rlt.isotope("Lu177").lambdaPhys;
            l = lp + [0.01 0.12 1.2];
            t = [0.25; 0.5; 1; 2; 4; 8; 24; 48; 96; 168];
            A = 2*exp(-l(1)*t) + 10*exp(-l(2)*t) + 60*exp(-l(3)*t);
            f = rlt.fitTAC(t, A, Model="triexp");
            tc.verifyEqual(f.tia, 2/l(1) + 10/l(2) + 60/l(3), 'RelTol', 5e-2);
            tc.verifyTrue(issorted(f.lambda));
        end

        function fitNeverImpliesNegativeActivityOrSuperPhysicalHalfLife(tc)
            iso = rlt.isotope("Lu177");
            t = [4; 24; 72; 168];
            A = [50; 80; 60; 30];                 % rising then falling
            for model = ["mono","uptake"]
                f = rlt.fitTAC(t, A, Model=model);
                tc.verifyGreaterThanOrEqual(f.lambda, iso.lambdaPhys - 1e-12);
                tc.verifyGreaterThanOrEqual(f.amplitude, 0);
                tc.verifyGreaterThanOrEqual(f.predict(0:10:300), 0);
            end
        end

        function tooFewPointsIsAnError(tc)
            tc.verifyError(@() rlt.fitTAC([24; 96; 168], [5; 3; 1], Model="biexp"), ...
                'rlt:tooFewPoints');
            tc.verifyError(@() rlt.fitTAC([1; 2; 3; 4; 5], [5; 4; 3; 2; 1], Model="triexp"), ...
                'rlt:tooFewPoints');
        end

        function badSeriesAreRejected(tc)
            tc.verifyError(@() rlt.fitTAC([24; 24; 96], [5; 4; 1]), 'rlt:badSeries');      % repeated time
            tc.verifyError(@() rlt.fitTAC([96; 24], [5; 4]), 'rlt:badSeries');             % not increasing
            tc.verifyError(@() rlt.fitTAC([24; 96], [0; 0]), 'rlt:badSeries');             % all zero
            tc.verifyError(@() rlt.fitTAC([24; 96], [5; 4; 1]), 'rlt:badSeries');          % length mismatch
            tc.verifyError(@() rlt.fitTAC([24; 96], [5; -1]), ?MException);                  % negative activity
            tc.verifyError(@() rlt.fitTAC([24; 96], [5; NaN]), ?MException);                 % NaN
        end

        function extrapolatedFractionMatchesClosedForm(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.01;
            t = [24; 96];
            f = rlt.fitTAC(t, 1000*exp(-lam*t), Model="mono");
            expected = 1 - (exp(-lam*24) - exp(-lam*96));
            tc.verifyEqual(f.extrapolatedFraction, expected, 'RelTol', 1e-3);
        end

        function qualityWarningsAreRaised(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.01;
            f = rlt.fitTAC([24; 96], 1000*exp(-lam*[24; 96]), Model="mono");
            tc.verifyTrue(any(contains(f.warnings, "exactly determined")));
            tc.verifyTrue(any(contains(f.warnings, "extrapolated")));
            g = rlt.fitTAC([24; 48; 96], [10; 20; 40], Model="mono");               % rising: no washout
            tc.verifyTrue(any(contains(g.warnings, "No biological clearance")));
        end

        % ---- model selection ----------------------------------------------

        function modelSelectionFallsBackToMonoWhenDataStarved(tc)
            iso = rlt.isotope("Lu177");
            t = [24; 96];
            A = 1000 * exp(-(iso.lambdaPhys + 0.01) * t);
            [best, cand, fits] = rlt.selectModel(t, A);
            tc.verifyEqual(best.model, "mono");
            tc.verifyEqual(height(cand), 1);
            tc.verifyEqual(numel(fits), 1);
        end

        function aiccIsNaNUntilEnoughPoints(tc)
            t = [24; 48; 96];
            A = [100; 70; 40];
            tc.verifyTrue(isnan(rlt.fitTAC(t, A, Model="mono").aicc));     % n-p-1 = 0
            t = [4; 24; 48; 96];
            A = [100; 90; 70; 40];
            tc.verifyTrue(isfinite(rlt.fitTAC(t, A, Model="mono").aicc));  % n-p-1 = 1
        end

        function selectionListsEveryEligibleModel(tc)
            lp = rlt.isotope("Lu177").lambdaPhys;
            t = [0.5; 1; 2; 4; 24; 48; 72];
            A = 5*exp(-(lp+0.01)*t) + 50*exp(-(lp+0.5)*t);
            [~, cand] = rlt.selectModel(t, A, Models=["mono","biexp","triexp"]);
            tc.verifyEqual(cand.Model, ["mono";"biexp";"triexp"]);
            tc.verifyEqual(sum(cand.Selected), 1);
        end

        % ---- estimators ----------------------------------------------------

        function hybridMatchesTruthForDenseMonoExponential(tc)
            iso = rlt.isotope("Lu177");
            lam = iso.lambdaPhys + 0.01;
            t = (0:2:1200)';
            A = 1000 * exp(-lam * t);
            tc.verifyEqual(rlt.tiaHybrid(t, A), 1000/lam, 'RelTol', 1e-3);
        end

        function hybridPhysicalTailNeverShorterThanEffectiveTail(tc)
            t = [24; 96];
            A = [100; 40];
            tc.verifyGreaterThanOrEqual( ...
                rlt.tiaHybrid(t, A, Tail="physical"), rlt.tiaHybrid(t, A, Tail="effective"));
        end

        function hybridDecayOriginRecoversBloodLikeCurves(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.3;
            t = (0.5:0.5:600)';
            A = 5 * exp(-lam * t);                 % starts high, first sample at 0.5 h
            truth = 5 / lam;
            decay = rlt.tiaHybrid(t, A, Origin="decay");
            zero  = rlt.tiaHybrid(t, A, Origin="zero");
            tc.verifyEqual(decay, truth, 'RelTol', 5e-3);
            tc.verifyLessThan(zero, decay);        % "zero" origin underestimates falling curves
        end

        function singlePointIsLeastSensitiveNearOneOverLambda(tc)
            lamTrue = 0.015;
            lamAssumed = 1.2 * lamTrue;
            truth = 100 / lamTrue;
            relErr = @(ts) abs(rlt.tiaSinglePoint(ts, 100*exp(-lamTrue*ts), lamAssumed) / truth - 1);
            tOpt = rlt.optimalSinglePointTime(lamTrue);
            tc.verifyEqual(tOpt, 1/lamTrue, 'RelTol', 1e-12);
            tc.verifyLessThan(relErr(tOpt), 0.03);
            tc.verifyLessThan(relErr(tOpt), relErr(24));
            tc.verifyLessThan(relErr(tOpt), relErr(168));
        end

        function estimateDispatcherAgreesWithDirectCalls(tc)
            t = [4; 24; 72; 168];
            A = [90; 70; 40; 15];
            tc.verifyEqual(rlt.estimate("mono", t, A), rlt.fitTAC(t, A).tia, 'RelTol', 1e-12);
            tc.verifyEqual(rlt.estimate("hybrid", t, A), rlt.tiaHybrid(t, A), 'RelTol', 1e-12);
            tc.verifyError(@() rlt.estimate("stp", t, A), 'rlt:missingLambdaPop');
        end

        % ---- decay correction ---------------------------------------------

        function decayCorrectionRoundTrips(tc)
            t = [4; 24; 96];
            A = [10; 7; 3];
            back = rlt.undoDecayCorrection(t, rlt.applyDecayCorrection(t, A, "Lu177"), "Lu177");
            tc.verifyEqual(back, A, 'RelTol', 1e-12);
            tc.verifyLessThan(rlt.undoDecayCorrection(96, 3, "Lu177"), 3);
        end

        % ---- dose ----------------------------------------------------------

        function absorbedDoseLocalDepositionConstant(tc)
            % 1 MBq*h = 3.6e9 decays; E = 0.134 MeV; 1 MeV = 1.602176634e-13 J; 1 kg
            tc.verifyEqual(rlt.absorbedDose(1, 1), 3.6e9*0.134*1.602176634e-13, 'RelTol', 1e-9);
            tc.verifyEqual(rlt.absorbedDose(2, 0.5), 4*rlt.absorbedDose(1, 1), 'RelTol', 1e-12);
            tc.verifyEqual(rlt.absorbedDose(10, 1, SValue=0.3), 3);
            tc.verifyError(@() rlt.absorbedDose(1, 1, AbsorbedFraction=1.5), 'rlt:badFraction');
        end

        % ---- virtual cohorts ----------------------------------------------

        function populationTruthMatchesNumericalIntegral(tc)
            for preset = ["kidney","tumour","blood"]
                pop = rlt.simulatePopulation(Preset=preset, N=5, Seed=3);
                t = (0:0.25:6000)';
                numeric = trapz(t, pop.eval(t), 2);
                tc.verifyEqual(numeric, pop.trueTIA, 'RelTol', 5e-3);
                tc.verifyGreaterThanOrEqual(pop.eval(t), -1e-12);
                tc.verifyGreaterThanOrEqual(pop.lambda(:), rlt.isotope(pop.isotope).lambdaPhys);
                tc.verifyTrue(all(diff(pop.lambda, 1, 2) >= 0, 'all'));
            end
        end

        function bloodPresetStartsAtItsAmplitude(tc)
            pop = rlt.simulatePopulation(Preset="blood", N=10, Seed=4);
            tc.verifyFalse(pop.hasUptake);
            tc.verifyEqual(pop.eval(0), pop.amp, 'RelTol', 1e-12);
            kid = rlt.simulatePopulation(Preset="kidney", N=10, Seed=4);
            tc.verifyEqual(kid.eval(0), zeros(10, 1), 'AbsTol', 1e-15);
        end

        function populationIsReproducible(tc)
            a = rlt.simulatePopulation(N=20, Seed=9);
            b = rlt.simulatePopulation(N=20, Seed=9);
            tc.verifyEqual(a.trueTIA, b.trueTIA);
            tc.verifyError(@() rlt.simulatePopulation(Preset="nope"), 'rlt:unknownPreset');
        end

        % ---- uncertainty and benchmarks -----------------------------------

        function propagateTIAIsReproducibleAndCentred(tc)
            t = [24; 96];
            A = [800; 450];
            a = rlt.propagateTIA(t, A, Method="mono", N=100, Seed=5);
            b = rlt.propagateTIA(t, A, Method="mono", N=100, Seed=5);
            tc.verifyEqual(a.tia, b.tia);
            tc.verifyEqual(a.median, a.central, 'RelTol', 0.15);
            tc.verifyGreaterThan(a.ci95(2), a.ci95(1));
        end

        function evaluateSchedulesReturnsOneRowPerScheduleAndMethod(tc)
            pop = rlt.simulatePopulation(N=20, Seed=2);
            R = rlt.evaluateSchedules({[24 96], [4 24 96]}, pop, Methods=["mono","hybrid","stp"]);
            tc.verifyEqual(height(R), 6);
            tc.verifyTrue(all(R.Within10Pct >= 0 & R.Within10Pct <= 100));
        end

        function evaluateSchedulesHandlesBloodPopulation(tc)
            pop = rlt.simulatePopulation(Preset="blood", N=15, Seed=2);
            R = rlt.evaluateSchedules({[1 24], [0.5 2 24 48]}, pop, Methods=["mono","hybrid"]);
            tc.verifyEqual(height(R), 4);
            tc.verifyEqual(R.NFailed, zeros(4, 1));
        end

    end
end
