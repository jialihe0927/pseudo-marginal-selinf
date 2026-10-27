# PMMH for selective inference: experiment code

This repository contains finalized experiment code and numerical results for
pseudo-marginal selective inference. Each experiment is self-contained. The
two applications include sampling-efficiency tables, global-null calibration,
directional-power figures, and sampling-validation analyses and coordinate ECDFs;
the Gaussian file-drawer experiment includes its
paper figure. The two power figures use the shared `plot_style.R`.

- [Gaussian file drawer](experiments/file_drawer/README_fd.md)
- [CV-Lasso](experiments/cv_lasso/README_lasso.md)
- [Changepoint](experiments/changepoint/README_cp.md)

The CV-Lasso and changepoint applications are divided into sampling efficiency,
global-null inference, non-null inference, and sampling-validation directories;
all four directories are organized as runnable public code. Sampling validation
supports plotting and offline diagnostics from supplied full samples, with
sampling available only by explicit request. The supplied RDS files contain the frozen
paper results. Simulation
entry points write to a user-specified new file and refuse to overwrite an
existing result. Run commands from the corresponding experiment directory.

## Data and Git LFS

Git and Git LFS are required to obtain the two large CV-Lasso inference files:

- `experiments/cv_lasso/inference_null/inference_null_lasso.rds`
  (34,423,726 bytes, about 32.83 MiB).
- `experiments/cv_lasso/inference_nonnull/inference_nonnull_lasso.rds`
  (140,486,247 bytes, about 133.98 MiB).

Only these two exact paths use Git LFS; other RDS files use ordinary Git.
After installing Git and Git LFS, clone and configure the clone locally:

```sh
git clone "<repository-url>" pmmh-selective-inference
cd pmmh-selective-inference
git lfs install --local
git lfs pull
```

Replace `<repository-url>` with the published repository URL. Keep the downloaded
RDS files at the paths above; the existing plotting and inference commands then
apply without conversion. A small text file whose first line is
`version https://git-lfs.github.com/spec/v1` is an LFS pointer, not the complete
RDS: run `git lfs pull` from the clone if the data were not downloaded.
Do not assume that GitHub's Download ZIP includes the full LFS data.
Remote LFS upload and clean-clone download have not yet been verified.

