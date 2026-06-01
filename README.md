# epi-preprocessing

MATLAB pipeline for preprocessing and reconstructing randomly undersampled 3D-EPI data acquired on GE and Siemens 3T scanners. Successor to [fmri-prep](https://github.com/rextlfung/fmri-prep). Sequence design lives in [rand3depi](https://github.com/rextlfung/rand3depi).

---

## Dependencies

| Library | Purpose |
|---|---|
| [Pulseq](https://github.com/pulseq/pulseq) | Siemens sequence definition (`mr.opts`) |
| [TOPPE](https://github.com/toppeMRI/toppe) | GE sequence definition & `ift3` |
| [BART](https://mrirecon.github.io/bart/) | ESPIRiT sensitivity maps, coil compression, L1-SENSE reconstruction |
| [Orchestra](https://github.com/rextlfung/orchestra) *(private)* | Reading GE ScanArchive `.h5` files |
| [hmriutils](https://github.com/rextlfung/hmriutils) | EPI ramp-sample gridding, odd/even ghost correction |

---

## Repository structure

```
epi-preprocessing/
├── config.m              # All user-editable paths and tunable parameters (edit this first)
├── run_preprocessing.m   # Batch driver — calls preprocess() for each sequence in cfg.seqnames
├── preprocess.m          # Stage 1 — raw data → zero-filled k-space volume (one sequence)
├── run_bart.m            # Stage 2 driver — L1-wavelet + TV reconstruction via BART pics
├── run_cg_sense.m        # Stage 2 driver — CG-SENSE reconstruction (no BART dependency)
├── recon_frames.m        # Stage 2 core — shared frame-loop, smaps loading, NIfTI write
├── cg_sense.m            # CG-SENSE solver (unregularized conjugate gradient)
├── set_seq_paths.m       # Helper — populate per-sequence cfg.seqname / cfg.fn.* fields
└── process_smaps.m       # Helper — mask, crop, and resize sensitivity maps
```

`params.m` is not stored in the repository — it is loaded at runtime from `cfg.fn.params`, which points to the per-acquisition sequence directory in the data folder.

---

## Pipeline overview

The pipeline runs in two sequential stages. Edit `config.m` to point at your data before running anything.

### Stage 1 — `run_preprocessing.m` / `preprocess.m`

`run_preprocessing.m` is the batch driver: it reads `cfg.seqnames` and calls `preprocess()` once per sequence. `preprocess.m` converts raw scanner output for one sequence into a zero-filled, coil-compressed k-space volume ready for iterative reconstruction.

```
noise scan  ──► whitening matrix
                     │
GRE data    ──► whiten ──► PCA compression matrix (Nvcoils virtual coils)
                     │                 │
                     ▼                 │
             ESPIRiT smaps ◄───────────┤
             (saved to disk)           │
                                       │
cal data    ──► whiten ──► compress ──►┤──► odd/even phase offsets (a)
                                       │    kx trajectory (kxo / kxe)
                                       │
EPI data    ──► whiten ──► compress ──► grid (NUFFT) ──► phase-correct
                                                  │
                                                  ▼
                                        scatter into zero-filled volume
                                        [Nx, Ny, Nz, Nvcoils, Nframes]
                                                  │
                                                  ▼
                                          write to .mat via matfile
                                          (frame-by-frame, bounded memory)
```

**Key memory note:** `ksp_gre` and `ksp_epi_raw` are never in memory simultaneously. The output is written frame-by-frame using `matfile` so the full time series is never allocated.

**Checkpoint/resume:** `last_completed_frame` is written to the output matfile after each frame. If `preprocess.m` is interrupted and restarted, it detects the checkpoint and skips already-completed frames.

**Nvcoils selection:** `Nvcoils` is chosen automatically from the eigenvalue spectrum of the whitened GRE sample covariance. Components are retained until `cfg.cc_energy_thresh` of total variance is explained, subject to a minimum of `2R` (ensuring enough virtual coils for SENSE reconstruction).

A quick sanity-check reconstruction (RSS or matched-filter SENSE) of the first 6 steady-state frames is displayed at the end via `interactive4D`.

---

### Stage 2 — `run_bart.m` / `run_cg_sense.m`

Both drivers read the zero-filled k-space produced by Stage 1 and delegate to the shared `recon_frames.m`, passing a per-frame reconstruction function handle. Output is one NIfTI per sequence with voxel sizes and TR embedded in the header.

```
zero-filled k-space  ──► (load all frames, then parfor)
                                  │
sensitivity maps     ──────────── ┤
                                  ▼
                   run_bart.m: BART pics -R W:7:0:λ_l1 -R T:7:0:λ_tv -i 100 -S
                   run_cg_sense.m: unregularized CG-SENSE
                                  │
                                  ▼
                          img [Nx, Ny, Nz, Nframes]
                                  │
                                  ▼
                          NIfTI write (with spatial metadata)
```

**`run_bart.m`** uses BART `pics` with L1-wavelet (`-R W:7:0:λ_l1`) and total-variation (`-R T:7:0:λ_tv`) regularization over spatial dims. The combined terms trigger ADMM. The `-S` flag is required because the randomized EPI trajectory has no ACS region; iteration count is 100 to ensure convergence at R = 6. Output filename: `<seqname>_recon_bart_l1_r<λ_l1>_tv_r<λ_tv>.nii`.

**`run_cg_sense.m`** uses the built-in `cg_sense.m` solver — a bare conjugate-gradient SENSE loop with no regularization. No BART dependency. Output filename: `<seqname>_recon_cgs_i<N>.nii`.

---

### Sensitivity map helper — `process_smaps.m`

Called from `preprocess.m` (Stage 1) and, as a fallback, from `recon_frames.m` (Stage 2). Applies a four-step post-processing pipeline to the raw maps returned by `makeSmaps`:

1. **Support mask** — threshold the last ESPIRiT eigenvalue map; zero out background voxels.
2. **z-crop** — trim the GRE volume symmetrically in z to match the EPI slab FOV. Assumes shared isocenter; mismatched isocenters will cause spatial misregistration.
3. **Interpolation** — `imresize3` to the EPI acquisition grid `[Nx, Ny, Nz]`.
4. **Normalization** — divide by the root-sum-of-squares across coils so `sum(|s_c|²) ≤ 1` everywhere, matching the ESPIRiT convention expected by BART `pics`.

---

## Sequence parameters

| Parameter | EPI | GRE (sensitivity maps) |
|---|---|---|
| Resolution | 2.4 × 2.4 × 2.4 mm | 2.0 × 2.0 × 2.0 mm |
| Matrix | 90 × 90 × 60 | 108 × 108 × 108 |
| FOV | 21.6 × 21.6 × 14.4 cm | 21.6 × 21.6 × 21.6 cm |
| TE | 30 ms | ~1/Δf_fat + 0.8 ms (in-phase) |
| TR (shot) | volumeTR / Nshots | 6 ms |
| Volume TR | 800 ms | — |
| Flip angle | Ernst (~13°) | Ernst (~4°) |
| Acceleration R | 6 | — |
| Echo train length | 75 | — |
| Virtual coils | auto (≥ 2R) | — |
| Regularization λ | 0.005 | — |
| B0 | 3 T | 3 T |

---

## Quick start

1. Clone the repo and add dependencies to your MATLAB path.
2. Edit `config.m` — set `cfg.datdir` and `cfg.seqnames` (one or more sequence names). All per-sequence paths are resolved automatically by `set_seq_paths.m`.
3. Run Stage 1:
   ```matlab
   run('run_preprocessing.m')
   ```
   For a single sequence set `cfg.seqnames = {'caipi_ts'}`. For batch, add more entries:
   ```matlab
   cfg.seqnames = {'caipi', 'caipi_ts', 'pd', 'pd_acs'};
   ```
   Set `cfg.interactive = false` to suppress blocking visualization windows during batch runs.
4. Inspect the sanity-check reconstruction. Adjust `cfg.delay` (k-space center offset) and `cfg.threshold_mask` if needed.
5. Run Stage 2 (choose one):
   ```matlab
   run('run_bart.m')      % L1-wavelet + TV SENSE via BART pics
   run('run_cg_sense.m')  % Unregularized CG-SENSE, no BART needed
   ```
   Both scripts process all sequences in `cfg.seqnames`.
6. The final reconstruction is saved as a NIfTI with embedded voxel sizes and TR:
   - BART: `<datdir>/recon/<seqname>_recon_bart_l1_r<λ>.nii`
   - CG-SENSE: `<datdir>/recon/<seqname>_recon_cgs_i<N>.nii`

---

## Configuration reference

All tunable parameters live in `config.m`. Key fields:

| Field | Default | Description |
|---|---|---|
| `cfg.seqnames` | `{'caipi_ts'}` | Cell array of sequence names to process; add entries for batch runs |
| `cfg.cc_energy_thresh` | 0.95 | Fraction of coil-data variance to retain; auto-selects Nvcoils from GRE eigenvalue spectrum with a 2R lower bound |
| `cfg.delay` | −1 | k-space center offset (samples); adjust if ghost artifacts appear |
| `cfg.SENSEmethod` | `'bart'` | `'bart'` (ESPIRiT) or `'pisco'` |
| `cfg.threshold_mask` | 1 | ESPIRiT eigenvalue threshold for support mask |
| `cfg.lamb_l1` | 0.005 | L1-wavelet regularization weight λ for BART `pics` (`-R W:7:0:λ`) |
| `cfg.lamb_tv` | 0.005 | Total-variation regularization weight for BART `pics` (`-R T:7:0:λ`) |
| `cfg.Nframes` | 30 | Maximum frames to reconstruct in Stage 2; defaults to all available frames in the matfile if the file contains fewer |
| `cfg.doSENSE` | `true` | `false` falls back to root-sum-of-squares |
| `cfg.interactive` | `true` | `false` suppresses blocking `interactive4D` windows and eigenvalue figure |
| `cfg.showEPIphaseDiff` | `true` | Plot odd/even phase difference during calibration |
| `cfg.useOrchestra` | `true` | Use Orchestra library to read ScanArchive `.h5` files |

---

## Related repositories

- Sequence design: [rextlfung/rand3depi](https://github.com/rextlfung/rand3depi)
- Previous pipeline: [rextlfung/fmri-prep](https://github.com/rextlfung/fmri-prep)
<!-- TODO: add paper/preprint link once available -->