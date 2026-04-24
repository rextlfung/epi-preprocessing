%{
cg_sense.m — CG-SENSE reconstruction of all temporal frames.

Reads the zero-filled k-space produced by main.m, loads (or recomputes)
sensitivity maps, then runs iterative L1-regularised SENSE reconstruction
via BART pics on every frame.  Output is written to a NIfTI file.

Dependencies: BART (https://mrirecon.github.io/bart/)
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');   % Loads cfg struct
run('./params.m');   % Loads sequence/system parameters

%% ── Derived filenames ────────────────────────────────────────────────────────
% Separate recon data directory (GRE and EPI may live elsewhere from storage)
datdir   = '/mnt/storage/rexfung/20260409tap/recon/';
fn_epi   = fullfile(datdir, 'caipi_ts_epi_zf.mat');
fn_gre   = fullfile(datdir, 'gre.mat');
fn_smaps = fullfile(datdir, sprintf('smaps_%s.mat', cfg.SENSEmethod));
fn_recon = sprintf('%s_recon_cgs_l1_r%.4f.nii', fn_epi(1:end-11), cfg.lamb);

%% ── BART reconstruction command ──────────────────────────────────────────────
% -l1 : L1 wavelet regularisation
% -r  : regularisation weight λ
% -e  : ESPIRiT model (uses sensitivity maps)
bart_cmd = sprintf('pics -l1 -r%f -e', cfg.lamb);
fprintf('BART command: %s\n', bart_cmd);

%% ── Load zero-filled k-space ─────────────────────────────────────────────────
fprintf('Loading EPI k-space from %s...\n', fn_epi);
try
    kdata = matfile(fn_epi);  % Use matfile to stream frames without loading all at once
catch ME
    error('cg_sense: Cannot open k-space file ''%s''.\n  %s', fn_epi, ME.message);
end

%% ── Load or compute sensitivity maps ─────────────────────────────────────────
if exist(fn_smaps, 'file')
    fprintf('Loading precomputed sensitivity maps from %s\n', fn_smaps);
    load(fn_smaps, 'smaps_raw', 'emaps');
else
    fprintf('Sensitivity maps not found. Estimating via %s...\n', cfg.SENSEmethod);
    try
        load(fn_gre, 'ksp_gre');
    catch ME
        error('cg_sense: Cannot load GRE data ''%s''.\n  %s', fn_gre, ME.message);
    end

    tic
        [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
    toc
    save(fn_smaps, 'smaps_raw', 'emaps', '-v7.3');
end

%% ── Process sensitivity maps ─────────────────────────────────────────────────
smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
    Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, cfg.Nvcoils, ...
    cfg.SENSEmethod, cfg.threshold_mask);

% Uncomment if x-direction alignment between GRE and EPI is needed:
% smaps = flip(smaps, 1);

%% ── Frame-by-frame reconstruction ───────────────────────────────────────────
Nframes = cfg.Nframes;

% Validate that the k-space file actually contains the requested frame count
kdata_size = size(kdata, 'ksp_epi_zf');  % Does not load data
if kdata_size(5) < Nframes
    error('cg_sense: cfg.Nframes (%d) exceeds frames in file (%d).', ...
        Nframes, kdata_size(5));
end

img = zeros(Nx, Ny, Nz, Nframes, 'single');

fprintf('Reconstructing %d frames...\n', Nframes);
tic
for frame = 1:Nframes
    fprintf('  Frame %d / %d\n', frame, Nframes);

    % Stream one frame from disk to avoid loading the full time series
    data = squeeze(kdata.ksp_epi_zf(:, :, :, :, frame));

    try
        img(:, :, :, frame) = bart(bart_cmd, data, smaps);
    catch ME
        warning('cg_sense: BART failed on frame %d — skipping.\n  %s', frame, ME.message);
    end
end
toc

%% ── Visualisation ────────────────────────────────────────────────────────────
interactive4D(abs(img));
return;   % Stop here during interactive use; NIfTI write below is manual

%% ── Write NIfTI ─────────────────────────────────────────────────────────────
fprintf('Writing reconstruction to %s\n', fn_recon);
niftiwrite(abs(img), fn_recon);
