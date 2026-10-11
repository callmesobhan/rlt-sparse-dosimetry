# Red-marrow dosimetry: methods, assumptions, pitfalls

Red marrow is usually the dose-limiting organ in radioligand therapy, and it is the hardest organ to do well: it cannot be imaged as a distinct organ, uptake is low, it is spread through the whole skeleton, and the quantity you can measure (blood, or a vertebral VOI) is only a proxy for the quantity you need (activity in the marrow cavities).

This document says what the toolbox computes, what each number assumes, and where the literature says the assumptions fail. Numbers in square brackets refer to the reference list at the end.

> **Status.** The equations and constants were checked against paper abstracts and secondary citations. The full text of the EANM bone-marrow guideline (Hindorf et al., EJNMMI 2010, DOI 10.1007/s00259-010-1422-4) was **not** available to the author of this toolbox. Before clinical use, compare every formula here with that guideline and with your department's procedure.

## 1. What is computed

Total marrow dose is a sum of contributions, each a TIA times an S-value [5, 8]:

```
D_RM = D_self  +  D_remainder-of-body  +  sum(organs)  +  sum(tumours)
```

| Term | How `rlt` gets the TIA | S-value |
|---|---|---|
| Marrow self-dose, blood-based | `rlt.marrow.selfDose(Source="blood")` | local deposition built in, or `SelfS` |
| Marrow self-dose, image-based | `rlt.marrow.selfDose(Source="marrow")` from a vertebral VOI | same |
| Remainder of body | `rlt.marrow.remainderOfBody`: whole-body TIA minus the other sources | **you supply** |
| Organs, tumours | `rlt.fitTAC` per region | **you supply** |

`rlt.marrow.totalDose` adds them and `rlt.marrow.uncertainty` propagates their errors. The toolbox ships **no S-values** (they depend on phantom and source model); read yours with `rlt.io.readSValues`, which insists on a recorded reference.

## 2. Blood-based self-dose

Following the AAPM/Sgouros blood model [1, 8]:

```
C_RM(t)  =  RMBLR · C_blood(t)          RMBLR = RMECFF / (1 − Hct)
```

With Hct 0.47 and RMECFF 0.19, RMBLR is 0.36 [1]; `rlt.marrow.rmblr` reproduces this. For a plasma basis the ratio is RMECFF alone (a red-marrow-to-plasma ratio of 0.19 is also used in the literature [10]).

The blood curve is fitted with a decreasing exponential (`mono`, `biexp`, `triexp`, chosen by AICc where the data allow), integrated analytically, and converted to dose. Under the default local-deposition model

```
D_self = RMBLR · TIAconc[MBq h/g] · 1000 · 3.6e9 · E_β · 1.602e-13 · φ
```

the marrow mass **cancels**: a heavier marrow holds more activity but dilutes the same energy. Mass only matters for the reported marrow TIA (MBq h) and for cross-dose terms from mass-scaled S-values.

## 3. Image-based self-dose

`Source="marrow"` treats the curve as the activity concentration in a vertebral VOI (for example L2–L4), multiplied by `Scale / RecoveryCoefficient`. Defaults are 1 and 1, so you must justify both:

- **RecoveryCoefficient** from phantom measurements at your VOI size and reconstruction settings.
- **Scale**: marrow concentration relative to the VOI concentration (the VOI also contains trabecular bone and fat). Methods differ here, and the choice of vertebrae alone changes dose, as does the presence of skeletal metastases [16].
- **LesionFree=true** is an assertion that you checked the VOI against lesion maps. Without it, a warning is recorded, because in a simulation study sites with bone lesions overestimated activity by a median of about 2.3× (range 1.06–10×) [12].

## 4. Which method can you trust?

| Situation | Blood-based | Image-based |
|---|---|---|
| No bone metastases, no specific marrow uptake | reasonable, a standard method [1, 6, 8] | usable |
| Specific marrow uptake (reported for ¹⁷⁷Lu-DOTATATE [13]) | **underestimates** | preferred |
| Bone metastases | **underestimates**; also a poor predictor of toxicity in one PSMA cohort [11] | needs lesion-free VOIs [12] or tumour/marrow localisation [11] |
| Few or no late scans | still possible | weak |

Evidence that the two families disagree:
- ¹⁷⁷Lu-DOTATATE, same cohort: blood-based median 15–28 mGy/GBq versus image-based 42–64 mGy/GBq, depending on the exponential model [14].
- ¹⁷⁷Lu-PSMA-617: image-based 25.4 versus blood-based 17.2 mGy/GBq, with the ratio rising over time, which suggests specific binding [15].
- ¹⁷⁷Lu-PSMA-617: a blood-surrogate (aorta) estimate was not concordant with image-based estimates and did not correlate with blood-count changes [11].
- Radioimmunotherapy: methods varied by a factor of 1.8–2.8 within the same patients [6]; image-based dose predicted toxicity better in a small lymphoma series [4].

**Recommendation encoded in the toolbox:** compute both when you can, report the discrepancy, and do not present a blood-only number as the marrow dose for patients with bone disease. `selfDose` always records the assumption in its output.

## 5. Pitfalls and how the toolbox responds

| Pitfall | Evidence | Toolbox response |
|---|---|---|
| Blood curve model changes the dose almost 2× | mono 15.3, bi 20.5, tri 28.3 mGy/GBq in the same patients [14] | `selfDose` fits every model the data support, lists each dose in `candidates`, reports `modelSpread`, warns above 25% |
| Too few samples to choose a model | AICc needs n ≥ p+2 (6 for biexp, 8 for triexp) | `selectModel` refuses to over-fit and falls back to the simplest model; `evaluateSchedules` with `Preset="blood"` benchmarks sampling plans |
| Marrow-to-blood ratio is not constant | RMBLR rose with time in rats and patients [2, 3] | `RatioFcn` option integrates a time-dependent ratio; the demo shows the sensitivity |
| Decay-corrected data fed to a fit | silent TIA overestimate | `DecayCorrected` column is mandatory on import; `rlt.undoDecayCorrection` |
| Percent entered as fraction (Hct 45) | classic data-entry error | `rlt.marrow.rmblr` and `readPatientInfo` reject it |
| Remainder of body negative | sources exceed whole body: units, calibration or model error | `remainderOfBody` errors by default instead of silently clipping |
| Cross-dose terms silently missing | | `rlt.io.lookupS` errors on a missing S-value; `totalDose` needs every source listed |
| Common calibration error ignored | shared by all image-derived TIAs | `uncertainty(..., SharedTIACV=)` |
| Different centres, same data, different dose | within ±30% even with one published method [8] | `provenance` and the report record every option, assumption and software version |
| Late-timepoint information | bi-exponential blood dose unchanged by samples after 24 h in one cohort, while a day-6 image added information [14] | benchmark your own schedules before dropping time points |

## 6. Assumptions register

| Item | Value in code | Source | Status |
|---|---|---|---|
| RMECFF | 0.19 (option) | Sgouros 1993 [1] | confirmed from abstract. Derived for activity that stays extracellular |
| RMBLR at Hct 0.47 | 0.36 | [1] | reproduced by `rmblr` |
| Hct | **no default**; user-supplied | | |
| Reference red-marrow mass | 1170 g (M), 900 g (F) | ICRP 70/89, via Pichardo 2007 | confirmed via secondary citation; check ICRP 89 table |
| Blood density | 1.06 g/mL | common reference value | **not verified here** |
| ¹⁷⁷Lu half-life, mean β energy | 6.647 d, 0.134 MeV | | **not verified here**; check ICRP 107 / NNDC (`rlt.isotope`) |
| β absorbed fraction in marrow | φ = 1 (option) | | upper bound for the β self-term; marrow geometry in trabecular bone is more complex |
| Photons, cross-dose | not in self-dose | | must be added with S-values |
| Uncertainty CVs | **no defaults** | | you must state them |

## 7. Not covered

Voxel-level Monte Carlo and spongiosa micro-structure models [9]; active-marrow localisation with a second tracer [11]; cortical-bone and bone-surface dose; α-emitters; dose–response modelling. Several are natural next steps once patient data exist (see the roadmap in `DESIGN.md`).

## References

1. Sgouros G. [Bone marrow dosimetry for radioimmunotherapy: theoretical considerations](https://consensus.app/papers/details/a6f87ae89a845428b34aa538843c6b4e/?utm_source=claude_desktop). J Nucl Med 1993.
2. Hindorf C et al. [Time dependence of the activity concentration ratio of red marrow to blood](https://consensus.app/papers/details/638f6389380e535e9710ada9e1773a84/?utm_source=claude_desktop). Cancer 2002.
3. Schwartz J et al. [Bone marrow dosimetry using 124I-PET](https://consensus.app/papers/details/daf315f1e0b756ea96d8dba9765f9625/?utm_source=claude_desktop). J Nucl Med 2012.
4. Ferrer L et al. [Three methods assessing red marrow dosimetry in lymphoma patients treated with radioimmunotherapy](https://consensus.app/papers/details/85a121fb963d5108a4d04a6ae5fc8ee4/?utm_source=claude_desktop). Cancer 2010.
5. Gosewisch A et al. [Patient-specific image-based bone marrow dosimetry in Lu-177-Octreotate and Lu-177-PSMA-617 therapy](https://consensus.app/papers/details/29204aa807015ba08db41243be881623/?utm_source=claude_desktop). EJNMMI Res 2018.
6. Hindorf C et al. [Evaluation of methods for red marrow dosimetry based on patients undergoing radioimmunotherapy](https://consensus.app/papers/details/4ea260dfc9af5a60b931d35e0438860e/?utm_source=claude_desktop). Acta Oncol 2005.
7. Hindorf C, Glatting G, Chiesa C, Lindén O, Flux G. EANM Dosimetry Committee guidelines for bone marrow and whole-body dosimetry. Eur J Nucl Med Mol Imaging 2010, DOI 10.1007/s00259-010-1422-4. *Not read in full.*
8. Wessels B et al. [Bone marrow dosimetry using blood-based models for radiolabeled antibody therapy: a multiinstitutional comparison](https://consensus.app/papers/details/87d9fdb141db5ebcb4a08d7e416cc8a3/?utm_source=claude_desktop). J Nucl Med 2004.
9. Blakkisrud J et al. [SPECT/CT image-derived absorbed dose to red marrow correlates with hematologic toxicity in [177Lu]Lu-DOTATATE patients](https://consensus.app/papers/details/ec5d83fc443c5414bb9c3ead90093745/?utm_source=claude_desktop). J Nucl Med 2024.
10. Makris N et al. [PET/CT-derived whole-body and bone marrow dosimetry of 89Zr-cetuximab](https://consensus.app/papers/details/7385cbde35045184853780800dc1edd3/?utm_source=claude_desktop). J Nucl Med 2015.
11. Peterson AB et al. [Comparison of imaging-based bone marrow dosimetry methodologies in [177Lu]Lu-PSMA-617 RLT](https://consensus.app/papers/details/032a7ca49e5157e2b15ec1de4d1c1668/?utm_source=claude_desktop). EJNMMI Phys 2025.
12. Rumiantcev M et al. [Feasibility of automated image-based red bone marrow dosimetry for [177Lu]Lu-PSMA](https://consensus.app/papers/details/b389db1556fb5108a976ba87ec2ba295/?utm_source=claude_desktop). Cancers 2025.
13. Hemmingsson J et al. [Specific uptake in the bone marrow causes high absorbed red marrow doses during [177Lu]Lu-DOTATATE treatment](https://consensus.app/papers/details/2b578ff67d2650c4bef65f43fce30dd1/?utm_source=claude_desktop). J Nucl Med 2023.
14. Vallot D et al. [Impact of different models based on blood samples and images for bone marrow dosimetry after 177Lu-labeled somatostatin-receptor therapy](https://consensus.app/papers/details/47889f7d80e555b8bc3dc68e3abbbe79/?utm_source=claude_desktop). EJNMMI Phys 2024.
15. Grob D et al. [Bone marrow dosimetry in low volume mHSPC patients receiving Lu-177-PSMA therapy using SPECT/CT](https://consensus.app/papers/details/538b27986fc05457bbf6e3e84028f67a/?utm_source=claude_desktop). EJNMMI Phys 2024.
16. Hagmarker L et al. [Bone marrow absorbed doses and correlations with hematologic response during 177Lu-DOTATATE treatments are influenced by image-based dosimetry method and presence of skeletal metastases](https://consensus.app/papers/details/b851e280273b5bb29852b72143b23ddc/?utm_source=claude_desktop). J Nucl Med 2019.
17. Pichardo JC et al. J Nucl Med 2007 (source of the ICRP 70/89 marrow masses quoted above; secondary citation).
