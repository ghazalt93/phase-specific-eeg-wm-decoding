# Order of execution

This file lists the intended execution order for the main article-1 workflow.

Helper functions should not be run directly. Raw EEG files, behavioral files, processed datasets, and large generated outputs are not included in this public repository.

## 1. Configuration

Edit local paths if needed:

```text
config/get_project_config.m
```

## 2. Optional raw-data processing

Run this only if raw subject folders are available locally:

```matlab
RUN00_BUILD_PROCESSED_FEATURE_DATASET_FROM_RAW
```

This script builds the processed feature dataset expected by the main analysis workflow. By default, the processed dataset is saved under:

```text
outputs/new_article_outputs/_wm_ml/dataset.mat
```

If the processed dataset has already been created, this step can be skipped.

## 3. Main manuscript workflow

Run the full main workflow with:

```matlab
RUN_ALL_REPRODUCE_MAIN_RESULTS
```

This wrapper checks for `dataset.mat` and then runs the following scripts in order:

1. `preprocessing_qc/RUN04_STEP100_subject_median_channel_roi_statistics.m`
2. `preprocessing_qc/RUN05_STEP104_loso_fisher_channel_roi_decoding_K1to22.m`
3. `preprocessing_qc/RUN06_STEP104_collect_loso_decoding_outputs.m`
4. `preprocessing_qc/RUN07_STEP106_permutation_highlighted_models.m`
5. `preprocessing_qc/RUN08_STEP107B_maxperm_primary_channel_binary_model.m`

### What each step does

- **RUN04**: subject-level channel and scalp-region feature statistics.
- **RUN05**: main leakage-controlled LOSO Fisher-score decoding grid.
- **RUN06**: collection of STEP104/RUN05 task-level outputs.
- **RUN07**: permutation tests for highlighted models.
- **RUN08**: max-permutation correction for the primary channel-level binary model family.

## 4. Synthetic demo

To test code execution without private data:

```matlab
RUN_DEMO_WITH_SYNTHETIC_DATA
```

The demo creates a small synthetic dataset and runs the early workflow steps on synthetic data. It is not intended to reproduce manuscript results.

## 5. Figure scripts

After main outputs exist, figure scripts may be run manually:

```matlab
manuscript_figures/RUN10_STEP215_make_decoding_presentation_figures
```

The feature-inventory/QC figure script is optional:

```matlab
manuscript_figures/RUN09_STEP28B_FEATURE_FAMILY_RECLASSIFY_AND_QC
```

## 6. Helper functions

The following folders contain functions required by the scripts above and should remain on the MATLAB path:

```text
feature_extraction/
helper_functions/
```

These functions should not be run as standalone scripts unless their headers explicitly say so.

## 7. Optional and archived files

The following folders are not part of the main reproducibility path:

```text
optional_qc/
optional_exploratory/
archive_legacy_statistics/
archive_legacy_local_checks/
archive_legacy/
archive_future_connectivity_article/
```

Connectivity-specific and exploratory scripts are retained only for traceability or future work and are not part of the present manuscript workflow.
