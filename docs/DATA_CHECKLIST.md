# Data checklist: what to collect so patient data drops straight in

The toolbox is built so that the first real patient needs no code changes. These are the things that, in practice, make a dosimetry dataset unusable after the fact. Align the protocol with your department's procedure and the EANM guidance; the sampling times below are examples from published work, not recommendations.

## 1. Per patient (`rlt.io.readPatientInfo`)

| Field | Why | Pitfall |
|---|---|---|
| Pseudonymous `PatientID` | links all files | no names or record numbers in IDs, file names or report labels |
| `Sex`, `Mass_kg`, `Height_cm` | reference marrow mass, S-value scaling | use the mass on the day of therapy |
| `Hct` on the day of therapy | blood-based ratio RMBLR = RMECFF/(1−Hct) | enter as a fraction (0.45), not percent |
| `BoneMetastases` yes/no/unknown, ideally burden from baseline PSMA PET or bone scan | decides whether blood-based or image-based marrow dose is trustworthy [see MARROW.md §4] | a blank here makes later interpretation impossible |
| Prior therapies (chemo, ²²³Ra, external beam, previous RLT cycles) | marrow reserve | |
| Platelets, neutrophils, Hb, WBC at baseline, weekly during therapy, at nadir, 3 and 6 months | the endpoint you will correlate dose against | record dates, not just values |

## 2. Per cycle

- Administered activity, calibration time and residual activity in the syringe/line, to the same clock as the scans.
- Start and end time of administration. Define **t = 0 once** (start of administration is the usual choice) and compute all `Time_h` from it.

## 3. Blood (region name `blood`)

- Sample actual clock times, not nominal ones.
- Whole blood or plasma; the counter, its calibration, and a cross-calibration against the SPECT/CT sensitivity (a shared calibration error is the largest correlated uncertainty, see `SharedTIACV`).
- Counts converted to activity concentration **at the time of counting**, with the counting time recorded. Enter values as measured at the sample time: set `DecayCorrected=0`, or `1` if (and only if) you already corrected to t = 0.
- Number of samples: the model you can fit is limited by it. AICc needs n ≥ parameters+2, so 6 samples to compare a bi-exponential, 8 for a tri-exponential. Published schedules include 30 min, 80 min, 24, 48, 72 h [Gosewisch 2018] and up to 9 samples [Vallot 2024]. In one ¹⁷⁷Lu-DOTATATE cohort a bi-exponential blood dose did not change with samples after 24 h [Vallot 2024]; check this with `examples/demo_marrow_blood_schedules.m` on your own cohort before relying on it.
- Include at least one early sample (within the first hour) and at least two late ones.

## 4. Imaging (regions such as `vertebrae_L2L4`, `wholebody`, `kidney_L`, `lesion_01`)

- Scan clock times; scanner and reconstruction settings (iterations, subsets, filter, scatter and attenuation correction, resolution recovery).
- **Calibration**: sensitivity factor from a phantom with known activity, plus **recovery coefficients** for the VOI sizes you use. These feed `RecoveryCoefficient`.
- **VOI definition**: which vertebrae, sphere diameter or automatic spongiosa segmentation, and a lesion map. Keep the lesion map, because `LesionFree=true` is only honest if you checked.
- One whole-body acquisition (planar or SPECT) at a late time point can pin down the remainder-of-body term [Gosewisch 2018]. Keep the count rates or whole-body activity, not just the derived TIA.
- Organ and tumour VOIs for every source you want in `totalDose`.

## 5. S-values (`rlt.io.readSValues`)

For each Source → red marrow pair you need: S, unit, isotope, the **reference** (software and version, phantom, mass), and ideally the phantom source mass. Never type S-values into a script without recording where they came from. A missing pair is an error by design.

## 6. File format rules (enforced by `rlt.io.readActivityTable`)

| Column | Rule |
|---|---|
| `PatientID, Cycle, Region` | one curve per combination |
| `Time_h` | hours after t = 0, finite, ≥ 0; no duplicates within a curve |
| `Value`, `Unit` | ≥ 0; Bq, kBq, MBq, (/mL or /g variants) |
| `DecayCorrected` | **required**, 0 or 1 |
| `Isotope` | optional; Lu177, Y90, I131 |

Templates with a synthetic example patient are in `data/templates/`. Run `rlt.selfCheck` once on every new machine.

## 7. Sanity checks before trusting a patient

1. `readActivityTable` reports no warnings you cannot explain.
2. For each fit, read `fit.warnings`: how much of the TIA is extrapolated, whether an exactly determined model was used.
3. `selfDose.modelSpread` and the `candidates` table: if the dose depends strongly on the model, say so in the report.
4. Whole-body TIA exceeds the sum of the source TIAs (`remainderOfBody` enforces this).
5. Blood-based and image-based marrow doses agree within your expectations. If not, look at bone disease and specific uptake before blaming the code.
6. Save the report (`rlt.report.marrow`) with the data file versions.
