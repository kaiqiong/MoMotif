# MoMotif  
**Modification of Motif Analysis at Single Base-Pair Resolution**

`MoMotif` is an R pipeline and reproducible workflow for identifying, refining, and validating **DNA binding motifs** from genome-wide assays such as ChIP-seq or ATAC-seq.  
It implements the computational framework introduced in  
**Lebeau *et al.*, Nucleic Acids Research, 2022**, providing modular tools for *differential binding*, *de novo motif discovery*, *motif extension*, *single-base resolution analysis*, and *biological validation*.

---

## 🧭 Overview

MoMotif (**Mo**dification of **Motif** analysis) enables precise, interpretable discovery and validation of transcription factor binding motifs through a modular pipeline.  
Each step can be executed independently or as part of a complete end-to-end workflow.

| Step | Description |
|------|--------------|
| **1. Differential binding** | Identify differentially bound genomic regions using [`csaw`](https://bioconductor.org/packages/csaw/). |
| **2. Motif discovery** | Apply **GADEM** for *de novo* motif discovery within enriched clusters. |
| **3. Motif extension** | Extend discovered motifs to capture flanking base-pair patterns and improve specificity. |
| **4. Matching status extraction** | Quantify motif–region matching across clusters or experimental conditions. |
| **5. Statistical validation** | Evaluate motif enrichment via chi-square or other statistical tests. |
| **6. Cross-context validation** | Validate motif activity across biological models (e.g., KI vs KIKI), additional datasets, or cancer cohorts. |
| **7. TAD comparison** | Relate motif distributions to topologically associating domain (TAD) boundaries for 3D genome context. |

---

## 📂 Repository Structure

```
vignette/
├── Step_1_differential_binding_using_csaw/
├── Step_2_GADEM_for_individual_clusters/
├── Step_3_Extend_motifs/
├── Step_4_extract_matching_status/
├── Step_5_chi_square_test_or_other/
├── Step_6_1_common_KI_KIKI/
├── Step_6_2_Other_validations/
├── Step_6_cancer_validation/
└── Step_7_TAD_comparison/
```



Each folder contains a self-contained R Markdown vignette documenting input, analysis, and output for the corresponding step.

---

📄 Citation

If you use MoMotif in your research, please cite:

Single base-pair resolution analysis of DNA binding motif with MoMotif reveals an oncogenic function of CTCF zinc-finger 1 mutation.
Benjamin Lebeau, Kaiqiong Zhao, Maika Jangal, Tiejun Zhao, Maria Guerra, Celia M T Greenwood, Michael Witcher.
Nucleic Acids Research, Volume 50, Issue 15, 26 August 2022, Pages 8441–8458.
https://doi.org/10.1093/nar/gkac658

