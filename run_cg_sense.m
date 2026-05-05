%{
run_cg_sense.m — Batch CG-SENSE reconstruction of all temporal frames.

Reads the zero-filled k-space produced by run_preprocessing.m, loads (or
recomputes) sensitivity maps, then runs iterative CG-SENSE reconstruction
via the built-in cg_sense() function on every frame.  Output is written to
a NIfTI file for each sequence in cfg.seqnames.  No external toolboxes
beyond standard MATLAB are required.

Sequences are processed in parallel across a parfor worker pool.
Set cfg.interactive = false (the default) when running in parallel —
graphics cannot be displayed from workers.

Dependencies: cg_sense.m (included in this repository)
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

parfor i = 1:numel(cfg.seqnames)
    recon_sequence(cfg, cfg.seqnames{i}, i, numel(cfg.seqnames));
end

fprintf('\nBatch complete.\n');

% ─────────────────────────────────────────────────────────────────────────────

function recon_sequence(cfg, seqname, idx, ntotal)

fprintf('\n[%d/%d] %s\n', idx, ntotal, seqname);

cfg = set_seq_paths(cfg, seqname);
run(cfg.fn.params);   % Loads MRI system + sequence parameters into workspace

%% ── Derived filenames ────────────────────────────────────────────────────────
% Separate recon data directory (GRE and EPI may live elsewhere from storage)
datdir   = strcat(cfg.datdir, 'recon/');
fn_epi   = fullfile(datdir, sprintf('%s_epi_zf.mat',              cfg.seqname));
fn_gre   = fullfile(datdir, 'gre.mat');
fn_smaps = fullfile(datdir, sprintf('smaps_%s.mat',               cfg.SENSEmethod));
fn_recon = fullfile(datdir, sprintf('%s_recon_cgs_i%d.nii',       cfg.seqname, cfg.num_iter));

%% ── Load zero-filled k-space ─────────────────────────────────────────────────
fprintf('Loading EPI k-space from %s...\n', fn_epi);
try
    kdata = matfile(fn_epi);  % Use matfile to stream frames without loading all at once
catch ME
    fprintf('ERROR: Cannot open k-space file ''%s''.\n  %s\nSkipping...\n', fn_epi, ME.message);
    return;
end

%% ── Load or compute sensitivity maps ─────────────────────────────────────────
if exist(fn_smaps, 'file')
    fprintf('Loading precomputed sensitivity maps from %s\n', fn_smaps);
    load(fn_smaps, 'smaps_raw', 'emaps', 'Nvcoils');
else
    fprintf('Sensitivity maps not found. Estimating via %s...\n', cfg.SENSEmethod);
    try
        load(fn_gre, 'ksp_gre');
    catch ME
        fprintf('ERROR: Cannot load GRE data ''%s''.\n  %s\nSkipping...\n', fn_gre, ME.message);
        return;
    end

    Nvcoils = size(ksp_gre, 4);  % infer from compressed GRE written by run_preprocessing.m
    tic
        [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
    toc
    save(fn_smaps, 'smaps_raw', 'emaps', 'Nvcoils', '-v7.3');
end

%% ── Process sensitivity maps ─────────────────────────────────────────────────
smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
    Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, ...
    cfg.SENSEmethod, cfg.threshold_mask);

% Uncomment if x-direction alignment between GRE and EPI is needed:
% smaps = flip(smaps, 1);

%% ── Frame-by-frame reconstruction ───────────────────────────────────────────
Nframes = cfg.Nframes;

% Validate that the k-space file actually contains the requested frame count
kdata_size = size(kdata, 'ksp_epi_zf');  % Does not load data
if kdata_size(5) < Nframes
    fprintf('ERROR: cfg.Nframes (%d) exceeds frames in file (%d). Skipping...\n', ...
        Nframes, kdata_size(5));
    return;
end

img = zeros(Nx, Ny, Nz, Nframes, 'single');

fprintf('Reconstructing %d frames with %d CG iterations...\n', Nframes, cfg.num_iter);
tic
for frame = 1:Nframes
    fprintf('  Frame %d / %d\n', frame, Nframes);

    % Stream one frame from disk to avoid loading the full time series
    data = squeeze(kdata.ksp_epi_zf(:, :, :, :, frame));

    try
        img(:, :, :, frame) = cg_sense(data, smaps, cfg.num_iter);
    catch ME
        warning('run_cg_sense: cg_sense failed on frame %d — skipping.\n  %s', frame, ME.message);
    end
end
toc

%% ── Visualisation ────────────────────────────────────────────────────────────
if cfg.interactive
    interactive4D(abs(img));
end

%% ── Write NIfTI ─────────────────────────────────────────────────────────────
fprintf('Writing reconstruction to %s\n', fn_recon);
niftiwrite(abs(img), fn_recon);

end
