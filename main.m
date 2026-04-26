%{
main.m — Load and preprocess randomised 3D EPI data.

Pipeline
────────
1.  Load configuration and sequence parameters.
2.  Load noise scan  →  whitening matrix.
3.  Load GRE data   →  whiten, PCA-compress, obtain compression matrix Vr.
4.  Estimate or load sensitivity maps (uses ksp_gre, then frees it).
5.  Load calibration data  →  odd/even phase offsets a, trajectory kxo/kxe.
6.  Load sampling log  →  schedules.
7.  Load EPI data (last, largest)  →  whiten, compress, grid, phase-correct,
    allocate zero-filled volume, save.
8.  Quick sanity-check reconstruction of a few frames.

Memory strategy: ksp_gre and ksp_epi_raw are never in memory simultaneously.
Each large intermediate is cleared the moment it is no longer needed.

Sequence design repository: https://github.com/rextlfung/rand3depi
%}

%% ── Dependencies ─────────────────────────────────────────────────────────────
run('./config.m');   % Loads cfg struct (paths, tunable parameters)
run('./params.m');   % Loads MRI system + sequence parameters into workspace

for p = cfg.addpaths
    addpath(p{1});
end

% Resolve kxoe filename now that Nx is available from params.m
cfg.fn.kxoe = fullfile(cfg.seqdir, sprintf('kxoe%d.mat', Nx));

%% ── Derived timing parameters ────────────────────────────────────────────────
NframesDiscard = round(discardDuration / volumeTR);  % Frames before steady state

%% STEP 1 — Noise scan  (small; needed to whiten every subsequent dataset)
fprintf('Loading noise scan...\n');
try
    ksp_noise = single(orc_read(cfg.fn.noise));
catch ME
    error('main: Failed to read noise file ''%s''.\n  %s', cfg.fn.noise, ME.message);
end

[Nfid, Ncoils, ~] = size(ksp_noise);
ksp_noise = permute(ksp_noise, [1 3 2]);
ksp_noise = reshape(ksp_noise, [], 1, 1, Ncoils);

%% STEP 2 — GRE data  →  cc_matrix (coil compression) + ksp_gre (for smaps)
fprintf('Loading GRE data...\n');
try
    ksp_gre_raw = single(orc_read(cfg.fn.gre));
catch ME
    error('main: Failed to read GRE file ''%s''.\n  %s', cfg.fn.gre, ME.message);
end

fprintf('  Max |Re(GRE)|: %g\n', max(real(ksp_gre_raw(:))));
fprintf('  Max |Im(GRE)|: %g\n', max(imag(ksp_gre_raw(:))));

% Reshape: discard the blipless calibration block prepended by the scanner
Nfid_gre = size(ksp_gre_raw, 1);
ksp_gre  = ksp_gre_raw(:, :, Nfid_gre+1:end);
ksp_gre  = reshape(ksp_gre, Nx_gre, Ncoils, Ny_gre, Nz_gre);
ksp_gre  = permute(ksp_gre, [1 3 4 2]);  % → [Nx_gre, Ny_gre, Nz_gre, Ncoils]
clear ksp_gre_raw;

% Whiten then PCA-compress to cfg.Nvcoils virtual coils.
% cc -A  : SVD over the full k-space volume (equivalent to the previous PCA approach)
% cc -M  : return the compression matrix rather than compressed data, so the
%          same matrix can be applied consistently to cal and EPI via ccapply.
ksp_gre   = bart('whiten -n', ksp_gre, ksp_noise);
cc_matrix = bart(sprintf('cc -p %d -A -M', cfg.Nvcoils), ksp_gre);
ksp_gre   = bart(sprintf('ccapply -p %d', cfg.Nvcoils), ksp_gre, cc_matrix);
save(strrep(cfg.fn.gre, '.h5', '.mat'), 'ksp_gre', '-v7.3');

%% STEP 3 — Sensitivity maps
% Done here so ksp_gre is freed BEFORE the large EPI data is loaded.
% Peak memory: ksp_noise + ksp_gre + smaps_raw  (no EPI yet).
fn_smaps = fullfile(cfg.datdir, sprintf('recon/smaps_%s.mat', cfg.SENSEmethod));

if cfg.doSENSE
    if exist(fn_smaps, 'file')
        fprintf('Loading precomputed sensitivity maps from %s\n', fn_smaps);
        load(fn_smaps, 'smaps_raw', 'emaps', 'smaps');
    else
        fprintf('Estimating sensitivity maps via %s...\n', cfg.SENSEmethod);
        tic
            [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
        toc
        smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
            Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, cfg.Nvcoils, ...
            cfg.SENSEmethod, cfg.threshold_mask);
        save(fn_smaps, 'smaps_raw', 'emaps', 'smaps', '-v7.3');
    end
    clear smaps_raw emaps;

    % Uncomment if x-direction alignment between GRE and EPI is needed:
    % smaps = flip(smaps, 1);
end

% plot smaps to check that they look reasonable
interactive4D(abs(smaps));

% ksp_gre is no longer needed — free it before loading EPI
clear ksp_gre;

%% STEP 4 — Calibration data  →  odd/even phase offsets a, trajectory kxo/kxe
% Small compared to EPI; process fully before loading EPI.
fprintf('Loading calibration data...\n');
try
    ksp_cal_raw = single(orc_read(cfg.fn.cal));
catch ME
    error('main: Failed to read calibration file ''%s''.\n  %s', cfg.fn.cal, ME.message);
end

assert(Nfid == size(ksp_cal_raw, 1), ...
    'main: Calibration Nfid (%d) does not match noise Nfid (%d).', ...
    size(ksp_cal_raw, 1), Nfid);
fprintf('  Max |Re(cal)|: %g,  Max |Im(cal)|: %g\n', ...
    max(real(ksp_cal_raw(:))), max(imag(ksp_cal_raw(:))));

% Whiten and coil-compress calibration data.
% Do NOT squeeze after whitening: keeping the singleton in dim 3 places coils
% in MATLAB dim 4 (= BART dim 3), which is where ccapply expects them.
ksp_cal = bart('whiten -n', ...
    reshape(permute(ksp_cal_raw, [1 3 2]), Nfid, [], 1, Ncoils), ksp_noise);
clear ksp_cal_raw;

ksp_cal = squeeze(bart(sprintf('ccapply -p %d', cfg.Nvcoils), ksp_cal, cc_matrix));
ksp_cal = permute(ksp_cal, [1 3 2]);  % → [Nfid, Nvcoils, N_cal]

% Load pre-computed k-space trajectory (cycles/cm) and apply gradient delay
try
    load(cfg.fn.kxoe, 'kxo', 'kxe');
catch ME
    error('main: Could not load kxoe file ''%s''.\n  %s', cfg.fn.kxoe, ME.message);
end
kxo = kxo / 100;
kxe = kxe / 100;

fprintf('Applying k-space center offset: %.2f samples\n', cfg.delay);
kxo = interp1(1:Nfid, kxo, (1:Nfid) - 0.5 - cfg.delay, 'linear', 'extrap');
kxe = interp1(1:Nfid, kxe, (1:Nfid) - 0.5 - cfg.delay, 'linear', 'extrap');

% Reshape calibration data and estimate odd/even phase offsets
% [a(1) = constant offset (rad), a(2) = linear term (rad/fov)]
ksp_cal  = reshape(permute(ksp_cal, [1 3 2]), Nfid, ETL, [], cfg.Nvcoils);
ETL_even = ETL - mod(ETL, 2);
oephase_data = hmriutils.epi.rampsampepi2cart( ...
    ksp_cal(:, 1:ETL_even, :, :), kxo, kxe, Nx, fov(1)*100, 'nufft');
oephase_data = ifftshift(ifft(fftshift(oephase_data), Nx, 1));
[a, ~]       = hmriutils.epi.getoephase( ...
    squeeze(mean(oephase_data, 3)), cfg.showEPIphaseDiff);
fprintf('  Constant phase offset: %.4f rad\n',     a(1));
fprintf('  Linear phase offset:   %.4f rad/fov\n', a(2));
clear ksp_cal oephase_data;

%% STEP 5 — Sampling log  (tiny; loaded here so it is ready for allocation)
try
    load(cfg.fn.samp_log, 'schedules');
catch ME
    error('main: Could not load sampling log ''%s''.\n  %s', cfg.fn.samp_log, ME.message);
end
[Nframes, Nshots, ETL, ~] = size(schedules);

%% STEP 6 — EPI data  (frame-by-frame to bound memory)
fprintf('Opening EPI archive: %s\n', cfg.fn.epi);
try
    epi_archive = GERecon('Archive.Load', cfg.fn.epi);
catch ME
    error('main: Failed to open EPI archive ''%s''.\n  %s', cfg.fn.epi, ME.message);
end

shots_per_frame       = ETL * Nshots;
total_shots_expected  = shots_per_frame * Nframes;
if epi_archive.FrameCount < total_shots_expected
    error('main: EPI archive has %d shots but %d expected (%d frames × %d shots/frame).', ...
        epi_archive.FrameCount, total_shots_expected, Nframes, shots_per_frame);
end

% Pre-allocate the output file so each frame can be written immediately,
% avoiding any need to hold the full time-series in memory.
fprintf('Pre-allocating output file: %s\n', cfg.fn.recon);
mf = matfile(cfg.fn.recon, 'Writable', true);
mf.ksp_epi_zf = complex(zeros(Nx, Ny, Nz, cfg.Nvcoils, Nframes, 'single'));

fprintf('Processing %d frames (%d shots/frame)...\n', Nframes, shots_per_frame);
tic
for frame = 1:Nframes
    fprintf('  Frame %d / %d\n', frame, Nframes);

    % ── Read ETL*Nshots consecutive readouts for this volumetric frame ─────
    ksp_frame_raw = zeros(Nfid, Ncoils, shots_per_frame, 'single');
    for s = 1:shots_per_frame
        shot = GERecon('Archive.Next', epi_archive);
        ksp_frame_raw(:, :, s) = single(shot.Data);
    end

    % ── Whiten ─────────────────────────────────────────────────────────────
    % Coils must be in BART dim 3 (MATLAB dim 4): [Nfid, shots, 1, Ncoils]
    ksp_frame = bart('whiten -n', ...
        reshape(permute(ksp_frame_raw, [1 3 2]), Nfid, shots_per_frame, 1, Ncoils), ...
        ksp_noise);
    clear ksp_frame_raw;

    % ── Coil-compress → [Nfid, cfg.Nvcoils, shots_per_frame] ──────────────
    ksp_frame = squeeze(bart(sprintf('ccapply -p %d', cfg.Nvcoils), ksp_frame, cc_matrix));
    ksp_frame = permute(ksp_frame, [1 3 2]);

    % ── Reshape → [Nfid, ETL, Nshots, Nvcoils] ────────────────────────────
    % Matches the layout the original produced via reshape+permute on the
    % full dataset before splitting across frames.
    ksp_frame = reshape(ksp_frame, Nfid, cfg.Nvcoils, ETL, Nshots);
    ksp_frame = permute(ksp_frame, [1 3 4 2]);

    % ── Grid along kx ──────────────────────────────────────────────────────
    ksp_frame_cart = hmriutils.epi.rampsampepi2cart( ...
        ksp_frame, kxo, kxe, Nx, fov(1)*100, 'nufft');

    % ── Odd/even phase correction ───────────────────────────────────────────
    ksp_frame_cart = hmriutils.epi.epiphasecorrect(ksp_frame_cart, a);

    % ── Scatter readouts into zero-filled volume ───────────────────────────
    ksp_frame_zf = zeros(Nx, Ny, Nz, cfg.Nvcoils, 'single');
    for shot_idx = 1:Nshots
        for echo = 1:ETL
            iy = schedules(frame, shot_idx, echo, 1);
            iz = schedules(frame, shot_idx, echo, 2);

            if any(ksp_frame_zf(:, iy, iz, :) ~= 0)
                warning('main: Overwriting frame %d, ky=%d, kz=%d. Check schedule.', ...
                    frame, iy, iz);
            end

            ksp_frame_zf(:, iy, iz, :) = ksp_frame_cart(:, echo, shot_idx, :);
        end
    end

    % ── Write this frame to disk immediately ───────────────────────────────
    mf.ksp_epi_zf(:, :, :, :, frame) = single(ksp_frame_zf);
end
toc

% ksp_noise and cc_matrix are no longer needed
clear ksp_noise cc_matrix;
fprintf('EPI processing complete. Output: %s\n', cfg.fn.recon);

%% STEP 7 — Quick sanity-check reconstruction
NtestFrames = 6;
% Stream only the test frames from disk; mf supports partial indexing
ksp_test = mf.ksp_epi_zf(:, :, :, :, NframesDiscard + (1:NtestFrames));

imgs_mc = zeros(Nx, Ny, Nz, cfg.Nvcoils, NtestFrames);
for frame = 1:NtestFrames
    imgs_mc(:, :, :, :, frame) = toppe.utils.ift3(ksp_test(:, :, :, :, frame));
end

% Coil combination
if cfg.doSENSE
    img_final = squeeze(sum(imgs_mc .* conj(smaps), 4));  % Matched-filter combination
else
    img_final = squeeze(sqrt(sum(abs(imgs_mc).^2, 4)));   % Root-sum-of-squares
end

interactive4D(abs(permute(img_final,   [2 3 1 4])));
interactive4D(angle(permute(img_final, [2 3 1 4])));