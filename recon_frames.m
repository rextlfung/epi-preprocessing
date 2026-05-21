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
        fprintf('ERROR: Cannot load GRE data ''%s''.\n  %s\nSkipping...\n', fn_gre, ME.message);
        return;
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

% Uncomment if x-direction alignment between GRE and EPI is needed:
% smaps = flip(smaps, 1);

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

%% Quick RSS reconstruction (zero-filled, no SENSE)
if cfg.interactive
    imgs_mc = zeros(Nx, Ny, Nz, Nvcoils, Nframes, 'single');
    for frame = 1:Nframes
        imgs_mc(:, :, :, :, frame) = toppe.utils.ift3(squeeze(ksp(:, :, :, :, frame)));
    end
    img_rss = squeeze(sqrt(sum(abs(imgs_mc).^2, 4)));
    interactive4D(abs(permute(img_rss, [2 3 1 4])));
    clear imgs_mc img_rss;
end

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
info = niftiinfo(fn_recon);
info.PixelDimensions = [[fov(1)/Nx, fov(2)/Ny, fov(3)/Nz] * 1e3, volumeTR];  % mm, s
info.SpaceUnits = 'Millimeter';
info.TimeUnits  = 'Second';
niftiwrite(abs(img), fn_recon, info);
end
