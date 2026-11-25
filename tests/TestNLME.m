classdef TestNLME < matlab.unittest.TestCase
    %TESTNLME Tests for Cycle-1 -> Cycle-2+ Empirical Bayes NLME optimizer and blood-proxy auditor.

    methods (TestClassSetup)
        function addRepoToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)

        function popPriorRecoversCohortDistribution(tc)
            stream = RandStream('mt19937ar', 'Seed', 11);
            iso = rlt.isotope("Lu177");
            t = [4; 24; 72; 168];
            N = 12;
            series = cell(N, 1);
            trueBio = 0.015 * exp(0.20 * randn(stream, N, 1));
            trueA0  = 100.0 * exp(0.25 * randn(stream, N, 1));
            for i = 1:N
                A = trueA0(i) * exp(-(iso.lambdaPhys + trueBio(i)) * t);
                series{i} = struct('t', t, 'A', A);
            end

            prior = rlt.nlme.fitPopPrior(series, Isotope="Lu177");
            tc.verifyEqual(prior.nPatients, N);
            tc.verifyEqual(prior.mu(1), mean(log(trueA0)), 'AbsTol', 0.05);
            tc.verifyEqual(prior.mu(2), mean(log(trueBio)), 'AbsTol', 0.08);
            tc.verifyGreaterThan(det(prior.omega), 0);
        end

        function cycle2BayesBeatsNaiveSTPWithCycle1Anchor(tc)
            stream = RandStream('mt19937ar', 'Seed', 19);
            iso = rlt.isotope("Lu177");
            tFull = [4; 24; 72; 168];
            N = 10;
            series = cell(N, 1);
            for i = 1:N
                lam = iso.lambdaPhys + 0.018 * exp(0.30 * randn(stream));
                a0  = 120 * exp(0.25 * randn(stream));
                series{i} = struct('t', tFull, 'A', a0 * exp(-lam * tFull));
            end
            prior = rlt.nlme.fitPopPrior(series, Isotope="Lu177");

            % Patient with slow biological clearance (outlier vs population): Cycle-1 anchor helps dramatically
            lamSlow = iso.lambdaPhys + 0.004;
            a0Slow  = 150.0;
            c1Fit = rlt.fitTAC(tFull, a0Slow * exp(-lamSlow * tFull), Isotope="Lu177", Model="mono");
            t2 = 24;  % early single time point where naive Hanscheid fails badly for slow clearers
            A2 = a0Slow * exp(-lamSlow * t2);
            trueTIA = a0Slow / lamSlow;

            resC1 = rlt.nlme.predictCycle2(t2, A2, prior, Cycle1Fit=c1Fit, ...
                MeasurementCV=0.05, DosePerTIA=1e-4, DoseLimitGy=0.5);
            errC1    = abs(resC1.tia / trueTIA - 1);
            errNaive = abs(resC1.tiaNaiveSTP / trueTIA - 1);
            tc.verifyLessThan(errC1, errNaive);
            tc.verifyLessThan(resC1.tiaCI95(1), resC1.tia);
            tc.verifyGreaterThan(resC1.tiaCI95(2), resC1.tia);
            tc.verifyTrue(resC1.probExceedLimit >= 0 && resC1.probExceedLimit <= 1);
        end

        function optimizeCycle2ScheduleFindsOptimalWindow(tc)
            stream = RandStream('mt19937ar', 'Seed', 23);
            iso = rlt.isotope("Lu177");
            tFull = [4; 24; 72; 168];
            series = cell(8, 1);
            for i = 1:8
                lam = iso.lambdaPhys + 0.012 * exp(0.25 * randn(stream));
                a0  = 100 * exp(0.20 * randn(stream));
                series{i} = struct('t', tFull, 'A', a0 * exp(-lam * tFull));
            end
            prior = rlt.nlme.fitPopPrior(series, Isotope="Lu177");

            out = rlt.nlme.optimizeCycle2Schedule(prior, {[24], [72], [120], [24 96]}, ...
                NVirtual=40, Seed=5);
            tc.verifyEqual(height(out.table), 12);
            tc.verifyTrue(strlength(out.bestSchedule) > 0);
            tc.verifyGreaterThanOrEqual(out.bestWithin10Pct, 0);
        end

        function auditBloodProxySeparatesValidFromSpecificMarrowUptake(tc)
            lp = rlt.isotope("Lu177").lambdaPhys;
            t = [2; 24; 48; 96; 168];
            rmblr = rlt.marrow.rmblr(0.45);
            cBlood = 2e-3 * exp(-(lp + 0.03) * t);

            % Case 1: Passive equilibrium (C_marrow == RMBLR * C_blood) -> VALID
            aValid = rlt.marrow.auditBloodProxy(t, cBlood, t, rmblr * cBlood, Hct=0.45);
            tc.verifyEqual(aValid.status, "VALID");
            tc.verifyEqual(aValid.doseRatio, 1.0, 'RelTol', 0.05);

            % Case 2: Specific marrow uptake / infiltration (slower clearance + 2.5x higher uptake) -> VIOLATED
            cInfilt = 2.2 * rmblr * 2e-3 * exp(-(lp + 0.005) * t);
            aBad = rlt.marrow.auditBloodProxy(t, cBlood, t, cInfilt, Hct=0.45);
            tc.verifyEqual(aBad.status, "VIOLATED");
            tc.verifyGreaterThan(aBad.doseRatio, 1.6);
            tc.verifyGreaterThan(aBad.rmblrSlopePerDay, 0.15);
            tc.verifyFalse(isempty(aBad.flags));
        end

    end
end
