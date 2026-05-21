%% config.m — Centralized configuration for the 3D EPI recon pipeline.
%
% Edit this file to update paths and tunable parameters without touching
% the pipeline scripts themselves.  Loaded via run('./config.m') at the
% top of run_preprocessing.m and cg_sense.m.  Per-sequence paths (cfg.fn.*, cfg.seqdir)
% are not set here — they are populated per-iteration by set_seq_paths.m.
%
% All settings are stored as fields of the struct `cfg` to avoid polluting
% the base workspace with ambiguous variable names.

%% ── Toolbox paths ────────────────────────────────────────────────────────────
cfg.addpaths = {
    '/home/rexfung/github/orchestra'   % Reading ScanArchives (GE private)
    '/home/rexfung/github/hmriutils'   % EPI odd/even ghost correction
    '/home/rexfung/github/toppe'       % TOPPE utilities (ift3, etc.)
};

%% ── Data directory & file names ──────────────────────────────────────────────
cfg.datdir   = '/StorageRAID/rexfung/20260409tap/';
cfg.seqnames = {'caipi', 'caipi_ts', 'pd'};  % Cell array of sequence names to process.
                              % Per-sequence paths are built by set_seq_paths.m.
cfg.fn.gre   = fullfile(cfg.datdir, 'scanarchives/gre.h5');

%% ── Coil parameters ─────────────────────────────────────────────────────────
% Nvcoils is auto-selected in preprocess.m from the whitened GRE eigenvalue spectrum.
% Components are kept until the cumulative explained variance reaches
% cfg.cc_energy_thresh. Lower bound: max(selected, 2*R) for SENSE feasibility.
cfg.cc_energy_thresh = 0.9;    % Retained portion of total coil-data variance.
                               % Lower → fewer virtual coils (faster, less SNR).
                               % Higher → more virtual coils (more SNR, slower).

%% ── EPI preprocessing ───────────────────────────────────────────────────────
cfg.delay            = -1;    % Estimated k-space center offset (samples).
                              % Negative values shift the echo earlier.
cfg.showEPIphaseDiff = true;  % Plot odd/even phase difference during calibration.

%% ── Sensitivity map estimation ───────────────────────────────────────────────
cfg.SENSEmethod    = 'bart';  % 'bart'  — uses BART ESPIRiT
                              % 'pisco' — uses PISCO (eigenvalue method)
cfg.threshold_mask = 1;       % Voxels whose last eigenvalue exceeds this
                              % threshold are zeroed out in the support mask.

%% ── CG-SENSE / PICS reconstruction ──────────────────────────────────────────
cfg.lamb_l1  = 0.005; % L1-wavelet regularisation weight (λ) for BART pics (-R W:7:0:λ).
                      % Larger → smoother, smaller → noisier but sharper.
cfg.lamb_tv  = 0.005; % Total-variation regularisation weight for BART pics (-R T:7:0:λ).
                      % Applied over spatial dims only (per-frame recon has no temporal dim).
cfg.num_iter = 100;   % Max CG iterations for cg_sense() (run_cg_sense.m).
                      % Matches the -i 100 flag used in run_bart.m.
cfg.Nframes  = Inf;   % Number of temporal frames to reconstruct. (Inf == all frames).

%% ── Pipeline options ─────────────────────────────────────────────────────────
cfg.useOrchestra = true;   % Use Orchestra library to read ScanArchive files.
cfg.doSENSE      = true;   % Estimate sensitivity maps and use SENSE combination.
                           % Set false to fall back to root-sum-of-squares.
cfg.useParfor    = false;  % Use parfor for parallel frame reconstruction.
                           % Set false to run serially (easier debugging).
cfg.interactive  = false;  % Set false to suppress blocking interactive4D calls
                           % and the eigenvalue figure (useful for batch runs).