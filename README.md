# Macrophage Reactivity Signature Defines Developmental, Sex, and Disease-Dependent States in the Lung.

## Summary
This repository contains notebooks for a multi-cohort analysis of lung macrophage reactivity using the SMaRT framework, with the C13 signature (note: lower C13 score = higher reactivity). We integrate microarray, bulk RNA-seq, and single-cell RNA-seq, from human and mouse lung to quantify how macrophage reactivity varies across development, sex, environmental exposure (smoking), and disease.

We found:
- Subset hierarchy: Interstitial macrophages (IM) exhibit higher reactivity than alveolar macrophages (AM) across diverse cohorts and technologies.
- Development & sex: Macrophage reactivity follows age-dependent trajectories and shows sex-linked differences consistent across platforms.
- Smoking exposure: Current smoking associates with reduced reactivity (higher C13 scores). 

## Requirements
- Python
  - [ScanPy](https://scanpy.readthedocs.io/en/stable/) (v1.9.1)
  - [Matplotlib](https://matplotlib.org/) (v3.5.3)
  - [Seaborn](https://seaborn.pydata.org/index.html) (v0.12.2)
  - [BoNE](https://github.com/sahoo00/BoNE) (Boolean Network Explorer)


## Citation
```
@article{lungDev2026,
  title  = {Macrophage Reactivity Signature Defines Developmental, Sex, and Disease-Dependent States in the Lung},
  author = {H M Zabir Haque, Betty Pham, Karen K Mestan, and Debashis Sahoo},
  note   = {Manuscript under review}
}
```
