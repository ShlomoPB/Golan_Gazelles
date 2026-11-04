# Golan mountain gazelle camera-trap project

This research repository contains the code, source data, derived data, model
objects, diagnostics, tables, and figures for the Golan Heights mountain
gazelle camera-trap study.

Reviewers should begin with [`REVIEWER_GUIDE.md`](REVIEWER_GUIDE.md), which
identifies the manuscript-critical inputs, scripts, and outputs and separates
them from exploratory and archival material.

## Repository contents

- Top-level `*.R` files and `R/`: preprocessing, spatial, temporal, demographic,
  sensitivity, visualization, and manuscript-support workflows.
- Core CSV and spreadsheet inputs: camera-round effort, site covariates,
  demographic count histories, classification corrections, and camera metadata.
- Analysis folders and saved outputs: model objects, diagnostics, predictions,
  figures, supplementary tables, and reproducibility reports.
- `student_data_all_species/`: date-corrected all-species exports prepared for
  supervised student use.

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
endangered species. Reuse or redistribution of coordinate-bearing files should
follow applicable permitting, data-ownership, and conservation-sensitivity
requirements.

No general reuse or redistribution license is granted by this repository.
