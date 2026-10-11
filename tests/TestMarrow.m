classdef TestMarrow < matlab.unittest.TestCase
    %TESTMARROW Bone-marrow dosimetry module tests.

    properties (Constant)
        DosePerMBqhPerKg = 3.6e9 * 0.134 * 1.602176634e-13   % Gy per (MBq h / kg), Lu-177 local deposition
    end

    methods (TestClassSetup)
        function addRepoToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)

        % ---- ratio and masses ---------------------------------------------

        function rmblrMatchesSgouros(tc)
            tc.verifyEqual(rlt.marrow.rmblr(0.47), 0.19/0.53, 'RelTol', 1e-12);
            tc.verifyEqual(rlt.marrow.rmblr(0.47), 0.36, 'AbsTol', 0.005);          % published value
            tc.verifyEqual(rlt.marrow.rmblr(Basis="plasma"), 0.19);
            tc.verifyEqual(rlt.marrow.rmblr(0.30, RMECFF=0.25), 0.25/0.70, 'RelTol', 1e-12);
        end

        function rmblrRejectsPercentAndImplausibleValues(tc)
            tc.verifyError(@() rlt.marrow.rmblr(45), 'rlt:badHct');
            tc.verifyError(@() rlt.marrow.rmblr(0), 'rlt:badHct');
            tc.verifyError(@() rlt.marrow.rmblr(NaN), 'rlt:badHct');
            tc.verifyError(@() rlt.marrow.rmblr(0.45, RMECFF=1.2), 'rlt:badRMECFF');
            [~, msg] = rlt.marrow.rmblr(0.70);
            tc.verifyGreaterThan(strlength(msg), 0);
            [~, msg] = rlt.marrow.rmblr(0.45);
            tc.verifyEqual(strlength(msg), 0);
        end

        function referenceMasses(tc)
            tc.verifyEqual(rlt.marrow.referenceMass("M"), 1170);
            tc.verifyEqual(rlt.marrow.referenceMass(" f "), 900);
            tc.verifyError(@() rlt.marrow.referenceMass("X"), 'rlt:badSex');
        end

        % ---- self-dose -----------------------------------------------------

        function selfDoseMatchesClosedForm(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.02;
            t = [2; 24; 48; 96];
            C0 = 2e-3;                                    % MBq/g at t = 0
            res = rlt.marrow.selfDose(t, C0*exp(-lam*t), ConcentrationUnit="MBq/g", ...
                Hct=0.45, Model="mono");
            tiaConc = C0 / lam;                           % MBq h / g blood
            expected = (0.19/0.55) * tiaConc * 1000 * tc.DosePerMBqhPerKg;
            tc.verifyEqual(res.doseGy, expected, 'RelTol', 1e-3);
            tc.verifyEqual(res.ratio, 0.19/0.55, 'RelTol', 1e-12);
        end

        function selfDoseDoesNotDependOnMarrowMass(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.02;
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-lam*t);
            a = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", MarrowMassG=800);
            b = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", MarrowMassG=1500);
            tc.verifyEqual(a.doseGy, b.doseGy, 'RelTol', 1e-12);
            tc.verifyEqual(b.marrowTIA / a.marrowTIA, 1500/800, 'RelTol', 1e-12);
            r = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", Sex="F");
            tc.verifyEqual(r.massG, 900);
            tc.verifyTrue(isnan(a.dosePerGBq));
        end

        function selfDoseConvertsMillilitresToGrams(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.02;
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-lam*t);
            perMl = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/mL", Hct=0.45, Model="mono");
            perG  = rlt.marrow.selfDose(t, C/1.06, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono");
            tc.verifyEqual(perMl.doseGy, perG.doseGy, 'RelTol', 1e-10);
        end

        function selfDoseReportsPerGBq(tc)
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-0.02*t);
            r = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", AdministeredMBq=7400);
            tc.verifyEqual(r.dosePerGBq, r.doseGy / 7.4, 'RelTol', 1e-12);
        end

        function selfDoseWithSValueMatchesLocalDepositionOfEqualS(tc)
            % An S-value equal to the local-deposition S for a 1000 g marrow must give the same dose.
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-0.02*t);
            S1000 = rlt.absorbedDose(1, 1.0);             % Gy per MBq h for 1 kg
            a = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono");
            b = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", ...
                SelfS=S1000, SelfSRefMassG=1000);
            tc.verifyEqual(b.doseGy, a.doseGy, 'RelTol', 1e-9);
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, SelfS=1e-5), ...
                'rlt:needRefMass');
        end

        function selfDoseInputGuards(tc)
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-0.02*t);
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g"), 'rlt:needHct');
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Ratio=0.3, ...
                RatioFcn=@(x) 0.3*ones(size(x))), 'rlt:badRatio');
            tc.verifyError(@() rlt.marrow.selfDose(t, C, Source="marrow"), 'rlt:needDensity');
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, ...
                Models=["mono","uptake"]), 'rlt:badModels');
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, ...
                AbsorbedFraction=1.2), 'rlt:badFraction');
            tc.verifyError(@() rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, ...
                RatioFcn=@(x) -ones(size(x))), 'rlt:badRatio');
        end

        function imageBasedSourceWarnsUnlessLesionFree(tc)
            t = [24; 72; 168];
            C = 5e-3 * exp(-0.015*t);
            r1 = rlt.marrow.selfDose(t, C, Source="marrow", ConcentrationUnit="MBq/g", Model="mono");
            tc.verifyTrue(any(contains(r1.warnings, "lesion-free")));
            r2 = rlt.marrow.selfDose(t, C, Source="marrow", ConcentrationUnit="MBq/g", Model="mono", ...
                LesionFree=true, RecoveryCoefficient=0.8);
            tc.verifyFalse(any(contains(r2.warnings, "lesion-free")));
            tc.verifyEqual(r2.ratio, 1/0.8, 'RelTol', 1e-12);
        end

        function timeDependentRatioIsIntegratedNumerically(tc)
            lam = rlt.isotope("Lu177").lambdaPhys + 0.03;
            t = [1; 4; 24; 48; 96];
            C = 1e-3 * exp(-lam*t);
            const = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Model="mono", Ratio=0.3);
            asFcn = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Model="mono", ...
                RatioFcn=@(x) 0.3*ones(size(x)));
            tc.verifyEqual(asFcn.doseGy, const.doseGy, 'RelTol', 1e-5);
            rising = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Model="mono", ...
                RatioFcn=@(x) 0.2 + 0.2*(1 - exp(-x/50)));
            flat = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Model="mono", Ratio=0.2);
            tc.verifyGreaterThan(rising.doseGy, flat.doseGy);
            tc.verifyGreaterThan(rising.ratioEffective, 0.2);
        end

        function candidatesExposeModelSpread(tc)
            lp = rlt.isotope("Lu177").lambdaPhys;
            t = [0.5; 1; 2; 4; 24; 48; 72];
            C = 2e-4 * (0.8*exp(-(lp+0.5)*t) + 0.17*exp(-(lp+0.08)*t) + 0.03*exp(-(lp+0.012)*t));
            r = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45);
            tc.verifyEqual(height(r.candidates), 3);
            tc.verifyEqual(sum(r.candidates.Selected), 1);
            tc.verifyGreaterThanOrEqual(r.modelSpread, 0);
            tc.verifyEqual(r.candidates.DoseGy(r.candidates.Selected), r.doseGy);
        end

        function elasticitiesMatchAnalyticValues(tc)
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-0.02*t);
            args = {"ConcentrationUnit", "MBq/g", "Model", "mono"};
            tc.verifyEqual(rlt.marrow.elasticity(t, C, "Hct", 0.45, args{:}), 0.45/0.55, 'AbsTol', 0.02);
            tc.verifyEqual(rlt.marrow.elasticity(t, C, "RMECFF", 0.19, "Hct", 0.45, args{:}), 1, 'AbsTol', 1e-6);
            tc.verifyEqual(rlt.marrow.elasticity(t, C, "AbsorbedFraction", 0.8, "Hct", 0.45, args{:}), 1, 'AbsTol', 1e-6);
        end

        % ---- remainder of body, S-values, totals --------------------------

        function remainderOfBodyArithmeticAndGuards(tc)
            [rob, info] = rlt.marrow.remainderOfBody(1000, [100 200 50]);
            tc.verifyEqual(rob, 650);
            tc.verifyEqual(info.fractionOfWholeBody, 0.35, 'AbsTol', 1e-12);
            tc.verifyError(@() rlt.marrow.remainderOfBody(100, [80 40]), 'rlt:negativeROB');
            [rob, info] = rlt.marrow.remainderOfBody(100, [80 40], Negative="clip");
            tc.verifyEqual(rob, 0);
            tc.verifyTrue(info.clipped);
        end

        function scaleSValueFollowsInverseMass(tc)
            tc.verifyEqual(rlt.marrow.scaleSValue(2e-5, 70, 140), 1e-5, 'RelTol', 1e-12);
            tc.verifyEqual(rlt.marrow.scaleSValue(2e-5, 70, 140, Exponent=0), 2e-5);
        end

        function totalDoseSumsContributions(tc)
            c = table(["Self";"Kidneys";"Remainder"], [100;50;400], [1e-3;2e-4;5e-5], ...
                'VariableNames', {'Source','TIA_MBqh','S_GyPerMBqh'});
            r = rlt.marrow.totalDose(c);
            tc.verifyEqual(r.totalGy, 0.1 + 0.01 + 0.02, 'RelTol', 1e-12);
            tc.verifyEqual(sum(r.table.SharePct), 100, 'AbsTol', 1e-9);
            dup = c;  dup.Source(2) = "Self";
            tc.verifyError(@() rlt.marrow.totalDose(dup), 'rlt:badContrib');
            neg = c;  neg.TIA_MBqh(1) = -1;
            tc.verifyError(@() rlt.marrow.totalDose(neg), 'rlt:badContrib');
        end

        function uncertaintyBehavesAtTheExtremes(tc)
            c = table(["Self";"Kidneys"], [100;50], [1e-3;2e-4], [0;0], [0;0], ...
                'VariableNames', {'Source','TIA_MBqh','S_GyPerMBqh','TIA_CV','S_CV'});
            u0 = rlt.marrow.uncertainty(c, N=200);
            tc.verifyEqual(u0.ci95, [0.11 0.11], 'RelTol', 1e-12);
            c.TIA_CV = [0.2; 0.2];  c.S_CV = [0.1; 0.1];
            u1 = rlt.marrow.uncertainty(c, N=4000, Seed=3);
            tc.verifyEqual(u1.median, 0.11, 'RelTol', 0.05);
            tc.verifyGreaterThan(u1.relHalfWidth, 0);
            tc.verifyEqual(sum(u1.perSource.VarianceSharePct), 100, 'AbsTol', 1e-6);
            u2 = rlt.marrow.uncertainty(c, N=4000, Seed=3, SharedTIACV=0.15);
            tc.verifyGreaterThan(u2.relHalfWidth, u1.relHalfWidth);
            tc.verifyError(@() rlt.marrow.uncertainty(c(:, 1:3)), 'rlt:missingUncertainty');
        end

        % ---- report --------------------------------------------------------

        function reportIsWrittenWithAssumptions(tc)
            t = [2; 24; 48; 96];
            C = 1e-3 * exp(-0.02*t);
            r = rlt.marrow.selfDose(t, C, ConcentrationUnit="MBq/g", Hct=0.45, Model="mono", ...
                AdministeredMBq=7400, Sex="M");
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(tmp.Folder, "report.md");
            rlt.report.marrow(file, r, Patient="TEST");
            txt = string(fileread(file));
            tc.verifyTrue(contains(txt, "Red-marrow dosimetry report"));
            tc.verifyTrue(contains(txt, "Assumptions"));
            tc.verifyTrue(contains(txt, "RMBLR"));
        end

    end
end
