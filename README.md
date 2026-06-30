# Phase-specific EEG decoding of visual working-memory conditions

This repository contains MATLAB code associated with the manuscript:

**Phase-specific channel and scalp-region EEG features decode visual working-memory conditions**

Manuscript ID: **ISCIENCE-D-26-08268**

## Overview

The code supports phase-specific EEG feature extraction, quality control, subject-level statistics, leakage-controlled feature selection, leave-one-subject-out decoding, permutation testing, and manuscript figure generation for a visual working-memory EEG study.

The task included color, orientation, and conjunction conditions across stimulus, maintenance, and retrieval phases. The final analyzed dataset included 22 participants, and the main decoding analyses used 22 leave-one-subject-out folds.

## Repository scope

The main workflow focuses on channel-level and scalp-region EEG features. Connectivity-based and exploratory analyses are not part of the main reproducibility path for this manuscript.

## Requirements

The code was developed with:

- MATLAB R2021a
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- EEGLAB, if EDF import or EEGLAB-compatible preprocessing is used
- BioSig, if EDF import relies on BioSig-supported routines

## Data availability

Raw EEG and behavioral files are not included in this repository because they contain human-participant data and are subject to institutional and ethical restrictions.

A processed, de-identified feature dataset may be made available by the corresponding author upon reasonable request and subject to institutional approval.

To reproduce the manuscript results, place the processed dataset at one of the expected locations under:

```text
outputs/new_article_outputs/
```

The default expected location is:

```text
outputs/new_article_outputs/_wm_ml/dataset.mat
```

## Configuration

Before running the analysis, check local paths in:

```text
config/get_project_config.m
```

The configuration file defines the project root, data root, output root, and active MATLAB paths.

Channel-location information used for plotting and channel-label checks is stored in `resources/channel_locations/ipm2.ced.tsv`.

## Main workflow

To run the full main workflow after placing `dataset.mat` in the expected location:

```matlab
RUN_ALL_REPRODUCE_MAIN_RESULTS
```

The main workflow runs:

1. `preprocessing_qc/RUN04_STEP100_subject_median_channel_roi_statistics.m`
2. `preprocessing_qc/RUN05_STEP104_loso_fisher_channel_roi_decoding_K1to22.m`
3. `preprocessing_qc/RUN06_STEP104_collect_loso_decoding_outputs.m`
4. `preprocessing_qc/RUN07_STEP106_permutation_highlighted_models.m`
5. `preprocessing_qc/RUN08_STEP107B_maxperm_primary_channel_binary_model.m`

## Synthetic demo

A small synthetic dataset can be generated to test code execution and path dependencies:

```matlab
RUN_DEMO_WITH_SYNTHETIC_DATA
```

The synthetic demo is only for testing code execution. It does not reproduce the manuscript results.

## Building the processed dataset from raw files

If raw data are available locally, the processed feature dataset can be generated with:

```matlab
RUN00_BUILD_PROCESSED_FEATURE_DATASET_FROM_RAW
```

Raw subject folders should be placed under the configured `data/Subjects` path, or the paths should be edited in the configuration/script.

## Outputs

Main outputs are written under:

```text
outputs/new_article_outputs/
```

Large generated outputs should not be committed to the repository.

## Optional and archived scripts

Folders named `optional_*` contain diagnostic or supplementary utilities. Folders named `archive_*` contain legacy or exploratory scripts retained for traceability. These folders are not part of the main manuscript reproducibility workflow.
