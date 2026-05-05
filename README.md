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
├── config.m           # All user-editable paths and tunable parameters (edit this first)
├── run_preprocessing.m             # Batch driver — calls preprocess() for each sequence in cfg.seqnames
├── preprocess.m       # Stage 1 — raw data → zero-filled k-space volume (one sequence)
├── cg_sense.m         # Stage 2 — CG-SENSE (BART pics) reconstruction → NIfTI
├── set_seq_paths.m    # Helper — populate per-sequence cfg.seqname / cfg.fn.* fields
└── process_smaps.m    # Helper — mask, crop, and resize sensitivity maps
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

**Nvcoils selection:** `Nvcoils` is chosen automatically from the eigenvalue spectrum of the whitened GRE sample covariance. Components are retained until `cfg.cc_energy_thresh` of total variance is explained, subject to a minimum of `2R` (ensuring enough virtual coils for SENSE reconstruction).

A quick sanity-check reconstruction (RSS or matched-filter SENSE) of the first 6 steady-state frames is displayed at the end via `interactive4D`.

---

### Stage 2 — `cg_sense.m`

Reconstructs every temporal frame from the zero-filled k-space produced by `run_preprocessing.m` using L1-regularised SENSE via BART `pics`.

```
zero-filled k-space  ──► (stream frame-by-frame)
                                  │
sensitivity maps     ──────────── ┤
                                  ▼
                          BART pics -l1 -r λ -i 100 -S
                                  │
                                  ▼
                          img [Nx, Ny, Nz, Nframes]
                                  │
                                  ▼
                          interactive4D  →  NIfTI write
```

The `-S` flag (strict SENSE) is used because the randomised EPI trajectory has no ACS region. Iteration count is set to 100 (vs. the BART default of 30) to ensure convergence at R = 6.

---

### Sensitivity map helper — `process_smaps.m`

Called from both `preprocess.m` and `cg_sense.m`. Applies a four-step post-processing pipeline to the raw maps returned by `makeSmaps`:

1. **Support mask** — threshold the last ESPIRiT eigenvalue map; zero out background voxels.
2. **z-crop** — trim the GRE volume symmetrically in z to match the EPI slab FOV. Assumes shared isocenter; mismatched isocenters will cause spatial misregistration.
3. **Interpolation** — `imresize3` to the EPI acquisition grid `[Nx, Ny, Nz]`.
4. **Normalisation** — divide by the root-sum-of-squares across coils so `sum(|s_c|²) ≤ 1` everywhere, matching the ESPIRiT convention expected by BART `pics`.

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
| Regularisation λ | 0.005 | — |
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
   Set `cfg.interactive = false` to suppress blocking visualisation windows during batch runs.
4. Inspect the sanity-check reconstruction. Adjust `cfg.delay` (k-space center offset) and `cfg.threshold_mask` if needed.
5. Run Stage 2:
   ```matlab
   run('cg_sense.m')
   ```
   `cg_sense.m` processes `cfg.seqnames{1}` by default; edit `cfg.seqnames{1}` or call `set_seq_paths` manually to target a different sequence.
6. The final reconstruction is saved as a NIfTI to `<datdir>/recon/<seqname>_recon_cgs_l1_r<λ>.nii`.

---

## Configuration reference

All tunable parameters live in `config.m`. Key fields:

| Field | Default | Description |
|---|---|---|
| `cfg.seqnames` | `{'caipi_ts'}` | Cell array of sequence names to process; add entries for batch runs |
| `cfg.cc_energy_thresh` | 0.9 | Fraction of coil-data variance to retain; auto-selects Nvcoils from GRE eigenvalue spectrum with a 2R lower bound |
| `cfg.delay` | −1 | k-space center offset (samples); adjust if ghost artifacts appear |
| `cfg.SENSEmethod` | `'bart'` | `'bart'` (ESPIRiT) or `'pisco'` |
| `cfg.threshold_mask` | 1 | ESPIRiT eigenvalue threshold for support mask |
| `cfg.lamb` | 0.005 | L1 regularisation weight λ for BART `pics` |
| `cfg.Nframes` | 30 | Frames to reconstruct in `cg_sense.m` |
| `cfg.doSENSE` | `true` | `false` falls back to root-sum-of-squares |
| `cfg.interactive` | `true` | `false` suppresses blocking `interactive4D` windows and eigenvalue figure |
| `cfg.showEPIphaseDiff` | `true` | Plot odd/even phase difference during calibration |
| `cfg.useOrchestra` | `true` | Use Orchestra library to read ScanArchive `.h5` files |

---

## Related repositories

- Sequence design: [rextlfung/rand3depi](https://github.com/rextlfung/rand3depi)
- Previous pipeline: [rextlfung/fmri-prep](https://github.com/rextlfung/fmri-prep)
<!-- TODO: add paper/preprint link once available -->