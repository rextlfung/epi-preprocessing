function [img, seq_params, runtime_s] = recon_frames(cfg, fn_recon, recon_fn)
%RECON_FRAMES  Load smaps + k-space, reconstruct all frames, write NIfTI.
%
%   recon_frames(cfg, fn_recon, recon_fn)
%   [img, seq_params, runtime_s] = recon_frames(cfg, fn_recon, recon_fn)
%
%   cfg        - Config struct with per-sequence paths set (from set_seq_paths)
%   fn_recon   - Output NIfTI file path; pass '' to skip NIfTI write
%   recon_fn   - Handle: @(data, smaps) -> [Nx, Ny, Nz] image for one frame
%
%   img        - Complex single [Nx, Ny, Nz, Nframes] reconstructed image
%   seq_params - Struct with geometry scalars: Nx, Ny, Nz, fov, volumeTR, Nvcoils, Nframes
%   runtime_s  - Wall-clock seconds for the parfor reconstruction loop

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
    error('recon_frames: cannot open k-space file ''%s'': %s', fn_epi, ME.message);
end

%% Load or compute sensitivity maps
if exist(fn_smaps, 'file')
    smaps_check = whos('-file', fn_smaps, 'smaps');
    if ~isempty(smaps_check)
        fprintf('Loading precomputed sensitivity maps from %s\n', fn_smaps);
        load(fn_smaps, 'smaps', 'Nvcoils');
    else
        fprintf('Loading raw smaps from %s (reprocessing)...\n', fn_smaps);
        load(fn_smaps, 'smaps_raw', 'emaps', 'Nvcoils');
        smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
            Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, ...
            cfg.SENSEmethod, cfg.threshold_mask);
        clear smaps_raw emaps;
    end
else
    fprintf('Sensitivity maps not found. Estimating via %s...\n', cfg.SENSEmethod);
    try
        load(fn_gre, 'ksp_gre');
    catch ME
        error('recon_frames: cannot load GRE data ''%s'': %s', fn_gre, ME.message);
    end
    Nvcoils = size(ksp_gre, 4);
    tic
        [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
    toc
    smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
        Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, ...
        cfg.SENSEmethod, cfg.threshold_mask);
    save(fn_smaps, 'smaps_raw', 'emaps', 'smaps', 'Nvcoils', '-v7.3');
end

% GRE and EPI both use a negative x pre-phaser (ArbEPI/GRE.m and
% ArbEPI/lib/make_prephasers.m), so no x-direction flip is needed.

%% Frame-by-frame reconstruction
kdata_size    = size(kdata, 'ksp_epi_zf');
Nframes_avail = kdata_size(5);
Nframes       = min(cfg.Nframes, Nframes_avail);
if Nframes < Nframes_avail
    fprintf('Reconstructing %d of %d available frames (cfg.Nframes cap).\n', Nframes, Nframes_avail);
end

img = zeros(Nx, Ny, Nz, Nframes, 'single');
Nworkers = 0;
if cfg.useParfor; Nworkers = Inf; end
ksp = kdata.ksp_epi_zf(:, :, :, :, 1:Nframes);

fprintf('Reconstructing %d frames...\n', Nframes);
t_start = tic;
parfor (frame = 1:Nframes, Nworkers)
    data = ksp(:, :, :, :, frame);
    try
        img(:, :, :, frame) = recon_fn(data, smaps);
    catch ME
        warning('recon_frames: reconstruction failed on frame %d — skipping.\n  %s', frame, ME.message);
    end
end
runtime_s = toc(t_start);
fprintf('Reconstruction done in %.1f s.\n', runtime_s);

if ~any(img(:))
    warning('recon_frames: all output frames are zero — recon_fn may have failed on every frame.');
end

seq_params = struct('Nx', Nx, 'Ny', Ny, 'Nz', Nz, ...
                    'fov', fov, 'volumeTR', volumeTR, ...
                    'Nvcoils', Nvcoils, 'Nframes', Nframes);

if cfg.interactive
    interactive4D(abs(img));
end

if ~isempty(fn_recon)
    fprintf('Writing reconstruction to %s\n', fn_recon);
    fn_tmp = [tempname() '.nii'];
    niftiwrite(abs(img), fn_tmp);
    info = niftiinfo(fn_tmp);
    delete(fn_tmp);
    info.PixelDimensions = [[fov(1)/Nx, fov(2)/Ny, fov(3)/Nz] * 1e3, volumeTR];  % mm, s
    info.SpaceUnits = 'Millimeter';
    info.TimeUnits  = 'Second';
    niftiwrite(abs(img), fn_recon, info);
end
end
