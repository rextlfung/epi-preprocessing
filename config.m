%% config.m — Centralized configuration for the 3D EPI recon pipeline.
%
% Edit this file to update paths and tunable parameters without touching
% the pipeline scripts themselves.  Loaded via run('./config.m') at the
% top of main.m and cg_sense.m.
%
% All settings are stored as fields of the struct `cfg` to avoid polluting
% the base workspace with ambiguous variable names.

%% ── Toolbox paths ────────────────────────────────────────────────────────────
cfg.addpaths = {
    '/home/rexfung/github/orchestra'   % Reading ScanArchives (GE private)
    '/home/rexfung/github/hmriutils'   % EPI odd/even ghost correction
};

%% ── Data directory & file names ──────────────────────────────────────────────
cfg.datdir  = '/mnt/storage/rexfung/20260409tap/recon/';
cfg.seqdir  = fullfile(cfg.datdir, '../seqs/pd/');

cfg.fn.gre      = fullfile(cfg.datdir, 'gre.h5');
cfg.fn.cal      = fullfile(cfg.datdir, 'pd_cal.h5');
cfg.fn.noise    = fullfile(cfg.datdir, 'pd_noise.h5');
cfg.fn.epi      = fullfile(cfg.datdir, 'pd_epi.h5');
cfg.fn.kxoe     = fullfile(cfg.seqdir, sprintf('kxoe%d.mat', 90)); % updated below
cfg.fn.samp_log = fullfile(cfg.seqdir, 'samp_locs.mat');
cfg.fn.recon    = fullfile(cfg.datdir, 'pd_epi_zf.mat');

% kxoe filename depends on Nx; set after params.m is loaded if needed
% cfg.fn.kxoe = fullfile(cfg.seqdir, sprintf('kxoe%d.mat', Nx));

%% ── Coil parameters ─────────────────────────────────────────────────────────
cfg.Nvcoils = 18;   % Virtual coils after PCA compression.
                    % Chosen based on visual inspection of the singular-value
                    % "knee". Decrease for faster recon, increase for more SNR.

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
cfg.lamb    = 0.01;  % L1 regularisation weight (λ) passed to BART pics.
                     % Larger → smoother, smaller → noisier but sharper.
cfg.Nframes = 18;    % Number of temporal frames to reconstruct in cg_sense.m.

%% ── Pipeline options ─────────────────────────────────────────────────────────
cfg.useOrchestra = true;   % Use Orchestra library to read ScanArchive files.
cfg.doSENSE      = true;   % Estimate sensitivity maps and use SENSE combination.
                           % Set false to fall back to root-sum-of-squares.
