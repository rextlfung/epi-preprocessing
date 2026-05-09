function recon_frames(cfg, fn_recon, recon_fn)
%RECON_FRAMES  Load smaps + k-space, reconstruct all frames, write NIfTI.
%
%   recon_frames(cfg, fn_recon, recon_fn)
%
%   cfg      - Config struct with per-sequence paths set (from set_seq_paths)
%   fn_recon - Output NIfTI file path
%   recon_fn - Handle: @(data, smaps) -> [Nx, Ny, Nz] image for one frame

run(cfg.fn.params);   % Loads Nx, Ny, Nz, fov, fov_gre, Nx_gre, ... into workspace

datdir   = strcat(cfg.datdir, 'recon/');
fn_epi   = fullfile(datdir, sprintf('%s_epi_zf.mat', cfg.seqname));
fn_gre   = fullfile(datdir, 'gre.mat');
fn_smaps = fullfile(datdir, sprintf('smaps_%s.mat', cfg.SENSEmethod));

%% Load zero-filled k-space
fprintf('Loading EPI k-space from %s...\n', fn_epi);
try
    kdata = matfile(fn_epi);
catch ME
    fprintf('ERROR: Cannot open k-space file ''%s''.\n  %s\nSkipping...\n', fn_epi, ME.message);
    return;
end

%% Load or compute sensitivity maps
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
    Nvcoils = size(ksp_gre, 4);
    tic
        [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
    toc
    save(fn_smaps, 'smaps_raw', 'emaps', 'Nvcoils', '-v7.3');
end

smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
    Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, ...
    cfg.SENSEmethod, cfg.threshold_mask);

% Uncomment if x-direction alignment between GRE and EPI is needed:
% smaps = flip(smaps, 1);

%% Frame-by-frame reconstruction
Nframes    = cfg.Nframes;
kdata_size = size(kdata, 'ksp_epi_zf');
if kdata_size(5) < Nframes
    fprintf('ERROR: cfg.Nframes (%d) exceeds frames in file (%d). Skipping...\n', ...
        Nframes, kdata_size(5));
    return;
end

img = zeros(Nx, Ny, Nz, Nframes, 'single');
Nworkers = 0;
if cfg.useParfor; Nworkers = Inf; end
ksp = kdata.ksp_epi_zf(:, :, :, :, 1:Nframes);
fprintf('Reconstructing %d frames...\n', Nframes);
tic
parfor (frame = 1:Nframes, Nworkers)
    data = squeeze(ksp(:, :, :, :, frame));
    try
        img(:, :, :, frame) = recon_fn(data, smaps);
    catch ME
        warning('recon_frames: reconstruction failed on frame %d — skipping.\n  %s', frame, ME.message);
    end
end
toc

if cfg.interactive
    interactive4D(abs(img));
end

fprintf('Writing reconstruction to %s\n', fn_recon);
niftiwrite(abs(img), fn_recon);
end
