# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the pipeline

All commands are run inside MATLAB (no shell build/test system).

**Stage 1 — preprocessing (raw k-space → zero-filled volume):**
```matlab
run('run_preprocessing.m')
```

**Stage 2 — reconstruction (choose one):**
```matlab
run('run_bart.m')        % L1-regularised SENSE via BART pics
run('run_cg_sense.m')   % Custom CG-SENSE (no BART dependency)
```

**Single sequence, ad-hoc:**
```matlab
run('./config.m');
cfg_seq = set_seq_paths(cfg, 'caipi_ts');
preprocess(cfg_seq);
```

**Batch run without blocking windows** — set `cfg.interactive = false` in `config.m` before running.

## Architecture

### Two-stage pipeline

```
config.m  →  run_preprocessing.m  →  preprocess.m       (Stage 1)
                                          ↓
                                  <seqname>_epi_zf.mat
                                          ↓
                              run_bart.m  OR  run_cg_sense.m  (Stage 2)
                                          ↓ (both delegate to recon_frames.m)
                                  <seqname>_recon_*.nii
```

Stage 1 is compute- and I/O-heavy (whitening, coil compression, NUFFT gridding). Stage 2 is iterative reconstruction; both stages support `parfor` parallelism at the frame level via `cfg.useParfor`.

`recon_frames.m` contains the shared Stage 2 logic: load/compute smaps, stream k-space, run the parfor loop, write NIfTI. `run_bart.m` and `run_cg_sense.m` each set up a reconstruction function handle (`@(data, smaps) bart(...)` or `@(data, smaps) cg_sense(...)`) and call `recon_frames`.

### cfg struct

`config.m` defines a single `cfg` struct that carries all parameters. It is never run as a function — it is executed with `run('./config.m')` to populate `cfg` in the caller's workspace. `set_seq_paths(cfg, seqname)` adds per-sequence fields (`cfg.seqname`, `cfg.seqdir`, `cfg.fn.*`) and returns the augmented struct.

`params.m` is **not in this repo** — it lives under `<datdir>/seqs/<seqname>/params.m` and is loaded at runtime via `run(cfg.fn.params)`. It dumps MRI system and sequence scalars (e.g., `Nx`, `Ny`, `Nz`, `ETL`, `R`, `fov`, `volumeTR`) directly into the caller's workspace, not into `cfg`.

### Memory strategy in preprocess.m

`ksp_gre` and `ksp_epi_raw` are never in memory simultaneously — `ksp_gre` is cleared before EPI data is loaded. The zero-filled output `ksp_epi_zf` is written frame-by-frame using `matfile` so the full time series is never allocated. Stage 2 scripts use `matfile` to stream frames during reconstruction for the same reason.

`preprocess.m` writes `last_completed_frame` to the output matfile after each frame. On restart, it detects this checkpoint and skips already-processed frames instead of starting over.

### Coil compression

`Nvcoils` is selected automatically from the eigenvalue spectrum of the whitened GRE covariance: components are kept until cumulative variance ≥ `cfg.cc_energy_thresh`, subject to a minimum of `2R`. The compression matrix (`cc_matrix`) is computed once from GRE and applied consistently to calibration and EPI data via BART `ccapply`.

Sensitivity maps are cached in `recon/smaps_<method>.mat`. `preprocess.m` saves the fully post-processed `smaps` variable (alongside `smaps_raw` and `emaps`). `recon_frames.m` loads `smaps` directly when it is present, skipping `process_smaps`; it falls back to reprocessing `smaps_raw`/`emaps` for legacy cache files that predate this convention. Smaps are recomputed from scratch only when the file is missing or `Nvcoils` has changed.

### Parallelism

Both Stage 2 scripts use `parfor (frame = 1:Nframes, Nworkers)`. `Nworkers` is set to `Inf` (full pool) when `cfg.useParfor = true`, or `0` (serial) when false. The full k-space array is loaded before the `parfor` loop so workers receive array slices rather than a `matfile` handle.

`cfg.Nframes` acts as an upper cap: `recon_frames.m` reconstructs `min(cfg.Nframes, actual_frames_in_file)` frames. NIfTI output includes voxel sizes (from `fov/N × 1000` mm) and `volumeTR` (seconds) in the header, sourced from `params.m`.

### Interactive / batch mode

`cfg.interactive = false` suppresses all blocking calls: `interactive4D`, the eigenvalue spectrum figure, and the odd/even phase difference plot (passed as `cfg.showEPIphaseDiff && cfg.interactive` to `getoephase`). Set this flag for any non-interactive or parallel run.

## External dependencies

| Library | Required for |
|---|---|
| [BART](https://mrirecon.github.io/bart/) | ESPIRiT smaps, coil compression, `run_bart.m` |
| [Orchestra](https://github.com/rextlfung/orchestra) *(private)* | Reading GE ScanArchive `.h5` files |
| [hmriutils](https://github.com/rextlfung/hmriutils) | EPI NUFFT gridding, odd/even phase correction |
| [TOPPE](https://github.com/toppeMRI/toppe) | GE sequence params (`ift3`) |
| [Pulseq](https://github.com/pulseq/pulseq) | Siemens sequence params |

Paths for Orchestra and hmriutils are set in `cfg.addpaths` inside `config.m`.

## Expected data layout

```
<datdir>/
├── scanarchives/
│   ├── gre.h5
│   ├── <seqname>_noise.h5
│   ├── <seqname>_cal.h5
│   └── <seqname>_epi.h5
├── seqs/<seqname>/
│   ├── params.m          # loaded at runtime; not in repo
│   ├── kxoe<Nx>.mat      # kx trajectory for odd/even readouts
│   └── samp_locs.mat     # random sampling schedule
└── recon/                # created by pipeline
    ├── gre.mat
    ├── smaps_<method>.mat
    ├── <seqname>_epi_zf.mat
    └── <seqname>_recon_*.nii
```
