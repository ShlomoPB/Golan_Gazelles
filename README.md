# Golan mountain gazelle camera-trap project

This private research repository contains the code, source data, derived data,
model objects, diagnostics, tables, figures, and submission-audit materials for
the Golan Heights mountain gazelle camera-trap study.

## Repository contents

- Top-level `*.R` files and `R/`: preprocessing, spatial, temporal, demographic,
  sensitivity, visualization, and manuscript-support workflows.
- Core CSV and spreadsheet inputs: camera-round effort, site covariates,
  demographic count histories, classification corrections, and camera metadata.
- Analysis folders and saved outputs: model objects, diagnostics, predictions,
  figures, supplementary tables, and reproducibility reports.
- `student_data_all_species/`: date-corrected all-species exports prepared for
  supervised student use.
- `submission_audit_output/`: manuscript and supplementary-information audit
  artifacts. These files do not replace the authoritative analysis inputs.

## Reproducibility

The repository intentionally retains saved model objects and derived outputs so
that reported results can be traced without automatically rerunning expensive
models. Session-information files accompany the major workflows where available.
Scripts should be reviewed for input paths and package requirements before use.

Large model and binary files are stored with [Git LFS](https://git-lfs.com/).
Install Git LFS before cloning, then use a normal clone:

```text
git lfs install
git clone <repository-url>
```

## Sensitive ecological data

This repository contains exact camera coordinates and fine-scale records for an
endangered species. Keep the repository **private**. Do not make it public,
redistribute the coordinate-bearing files, or publish fine-scale locations
without a formal sensitivity, permitting, and data-ownership review. Public
release should use a separately prepared, spatially de-identified data package.

No general reuse or redistribution license is granted by this repository.
