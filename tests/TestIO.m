classdef TestIO < matlab.unittest.TestCase
    %TESTIO Data import and validation tests.

    methods (TestClassSetup)
        function addRepoToPath(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Access = private)
        function folder = tempFolder(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            folder = tmp.Folder;
        end

        function T = goodTable(~)
            T = table(["P1";"P1";"P1";"P1";"P1"], [1;1;1;1;1], ...
                ["blood";"blood";"blood";"wholebody";"wholebody"], ...
                [4;24;96;24;72], [10;5;1;300;150], ["kBq/mL";"kBq/mL";"kBq/mL";"MBq";"MBq"], [1;1;1;0;0], ...
                'VariableNames', {'PatientID','Cycle','Region','Time_h','Value','Unit','DecayCorrected'});
        end
    end

    methods (Test)

        function activityTableRoundTrip(tc)
            file = fullfile(tc.tempFolder(), "a.csv");
            writetable(tc.goodTable(), file);
            [T, issues] = rlt.io.readActivityTable(file);
            tc.verifyEqual(height(issues), 0);
            lp = rlt.isotope("Lu177").lambdaPhys;

            s = rlt.io.getSeries(T, "P1", 1, "blood");
            tc.verifyEqual(s.t, [4; 24; 96]);
            tc.verifyEqual(s.unit, "MBq/mL");
            tc.verifyEqual(s.A, [10; 5; 1] * 1e-3 .* exp(-lp*[4; 24; 96]), 'RelTol', 1e-12);   % kBq->MBq, decay undone

            w = rlt.io.getSeries(T, "P1", 1, "wholebody");
            tc.verifyEqual(w.A, [300; 150]);                % not decay-corrected: left as is
            tc.verifyEqual(w.unit, "MBq");
        end

        function unsortedRowsAreSorted(tc)
            T = tc.goodTable();
            T = T([3 1 2 4 5], :);
            file = fullfile(tc.tempFolder(), "a.csv");
            writetable(T, file);
            tc.verifyWarning(@() rlt.io.readActivityTable(file), 'rlt:tableWarning');
            w = warning('off', 'rlt:tableWarning');
            restore = onCleanup(@() warning(w));
            S = rlt.io.getSeries(rlt.io.readActivityTable(file), "P1", 1, "blood");
            tc.verifyTrue(issorted(S.t));
        end

        function missingDecayCorrectionFlagIsAnError(tc)
            T = removevars(tc.goodTable(), 'DecayCorrected');
            file = fullfile(tc.tempFolder(), "a.csv");
            writetable(T, file);
            tc.verifyError(@() rlt.io.readActivityTable(file), 'rlt:invalidTable');
            iss = rlt.io.validateActivityTable(T);
            tc.verifyTrue(any(contains(iss.Message, "DecayCorrected")));
        end

        function allProblemsAreReportedTogether(tc)
            T = tc.goodTable();
            T.Unit(1) = "mCi";                % unknown unit
            T.Value(2) = -3;                  % negative
            T.Time_h(3) = T.Time_h(2);        % duplicate time in the group
            T.DecayCorrected(4) = 7;          % invalid flag (also differs inside its group)
            iss = rlt.io.validateActivityTable(T);
            err = iss(iss.Severity == "error", :);
            tc.verifyGreaterThanOrEqual(height(err), 4);
            file = fullfile(tc.tempFolder(), "a.csv");
            writetable(T, file);
            try
                rlt.io.readActivityTable(file);
                tc.verifyFail("expected an error");
            catch ME
                tc.verifyEqual(string(ME.identifier), "rlt:invalidTable");
                tc.verifySubstring(string(ME.message), "mCi");
                tc.verifySubstring(string(ME.message), "Duplicate");
            end
        end

        function mixedUnitsAndDecayFlagsInOneCurveAreErrors(tc)
            T = tc.goodTable();
            T.Unit(2) = "MBq";
            iss = rlt.io.validateActivityTable(T);
            tc.verifyTrue(any(contains(iss.Message, "mixes units")));
            T = tc.goodTable();
            T.DecayCorrected(2) = 0;
            iss = rlt.io.validateActivityTable(T);
            tc.verifyTrue(any(contains(iss.Message, "mixes decay-corrected")));
        end

        function timesInDaysTriggerAWarning(tc)
            T = tc.goodTable();
            T.Time_h(1:3) = [0.1; 0.5; 1];     % a curve whose times all sit below 3 h
            iss = rlt.io.validateActivityTable(T);
            tc.verifyTrue(any(contains(iss.Message, "Are times in hours")));
        end

        function getSeriesListsAvailableRegions(tc)
            file = fullfile(tc.tempFolder(), "a.csv");
            writetable(tc.goodTable(), file);
            T = rlt.io.readActivityTable(file);
            try
                rlt.io.getSeries(T, "P1", 1, "kidney");
                tc.verifyFail("expected an error");
            catch ME
                tc.verifyEqual(string(ME.identifier), "rlt:seriesNotFound");
                tc.verifySubstring(string(ME.message), "wholebody");
            end
        end

        function patientInfoValidation(tc)
            folder = tc.tempFolder();
            P = table(["P1";"P2"], ["M";"f"], [0.45;0.40], [80;62], ["yes";"0"], ...
                'VariableNames', {'PatientID','Sex','Hct','Mass_kg','BoneMetastases'});
            file = fullfile(folder, "p.csv");
            writetable(P, file);
            R = rlt.io.readPatientInfo(file);
            tc.verifyEqual(R.Sex, ["M";"F"]);
            tc.verifyEqual(R.BoneMetastases, ["yes";"no"]);

            P.Hct(1) = 45;                                 % percent instead of fraction
            writetable(P, file);
            tc.verifyError(@() rlt.io.readPatientInfo(file), 'rlt:invalidTable');

            P.Hct(1) = 0.45;  P.PatientID(2) = "P1";       % duplicate id
            writetable(P, file);
            tc.verifyError(@() rlt.io.readPatientInfo(file), 'rlt:invalidTable');
        end

        function sValueTableConvertsUnitsAndDemandsProvenance(tc)
            folder = tc.tempFolder();
            S = table(["RedMarrow";"RedMarrow"], ["Kidneys";"RemainderOfBody"], ["Lu177";"lu177"], ...
                [1.0;2.0], ["mGy/(MBq s)";"Gy/(MBq h)"], ["Test phantom v1";"Test phantom v1"], ...
                'VariableNames', {'Target','Source','Isotope','S','Unit','Reference'});
            file = fullfile(folder, "s.csv");
            writetable(S, file);
            R = rlt.io.readSValues(file);
            tc.verifyEqual(R.S_GyPerMBqh, [3.6; 2.0], 'RelTol', 1e-12);
            tc.verifyEqual(rlt.io.lookupS(R, "RedMarrow", "Kidneys", "Lu177"), 3.6, 'RelTol', 1e-12);
            tc.verifyError(@() rlt.io.lookupS(R, "RedMarrow", "Liver", "Lu177"), 'rlt:missingSValue');

            S.Reference(1) = "";
            writetable(S, file);
            tc.verifyError(@() rlt.io.readSValues(file), 'rlt:invalidTable');

            S.Reference(1) = "PLACEHOLDER - fill in";
            writetable(S, file);
            tc.verifyWarning(@() rlt.io.readSValues(file), 'rlt:placeholderSValue');
        end

        function missingFilesAreReported(tc)
            tc.verifyError(@() rlt.io.readActivityTable("does_not_exist.csv"), 'rlt:fileNotFound');
            tc.verifyError(@() rlt.io.readPatientInfo("does_not_exist.csv"), 'rlt:fileNotFound');
            tc.verifyError(@() rlt.io.readSValues("does_not_exist.csv"), 'rlt:fileNotFound');
        end

        function shippedTemplatesAreValid(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tdir = fullfile(root, "data", "templates");
            [T, issues] = rlt.io.readActivityTable(fullfile(tdir, "activity_table_template.csv"));
            tc.verifyEqual(sum(issues.Severity == "error"), 0);
            tc.verifyGreaterThan(height(T), 0);
            P = rlt.io.readPatientInfo(fullfile(tdir, "patient_info_template.csv"));
            tc.verifyGreaterThan(height(P), 0);
            tc.verifyWarning(@() rlt.io.readSValues(fullfile(tdir, "svalues_template.csv")), ...
                'rlt:placeholderSValue');
        end

    end
end
