function preprocess(cfg)
%{
preprocess — Load and preprocess randomised 3D EPI data for one sequence.

Pipeline
────────
1.  Load sequence parameters.
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

Called from run_preprocessing.m (batch) or directly:
    preprocess(set_seq_paths(cfg, 'caipi_ts'))

Sequence design repository: https://github.com/rextlfung/rand3depi
%}

%% ── Dependencies ─────────────────────────────────────────────────────────────
run(cfg.fn.params);  % Loads MRI system + sequence parameters into workspace

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
    error('preprocess: Failed to read noise file ''%s''.\n  %s', cfg.fn.noise, ME.message);
end

[Nfid, Ncoils, ~] = size(ksp_noise);
ksp_noise = permute(ksp_noise, [1 3 2]);
ksp_noise = reshape(ksp_noise, [], 1, 1, Ncoils);

%% STEP 2 — GRE data  →  cc_matrix (coil compression) + ksp_gre (for smaps)
fprintf('Loading GRE data...\n');
try
    ksp_gre_raw = single(orc_read(cfg.fn.gre));
catch ME
    error('preprocess: Failed to read GRE file ''%s''.\n  %s', cfg.fn.gre, ME.message);
end

fprintf('  Max |Re(GRE)|: %g\n', max(real(ksp_gre_raw(:))));
fprintf('  Max |Im(GRE)|: %g\n', max(imag(ksp_gre_raw(:))));

% Reshape: discard the blipless calibration block prepended by the scanner
Nfid_gre = size(ksp_gre_raw, 1);
ksp_gre  = ksp_gre_raw(:, :, Ny_gre+1:end);
ksp_gre  = reshape(ksp_gre, Nx_gre, Ncoils, Ny_gre, Nz_gre);
ksp_gre  = permute(ksp_gre, [1 3 4 2]);  % → [Nx_gre, Ny_gre, Nz_gre, Ncoils]
clear ksp_gre_raw;

% Whiten, auto-select Nvcoils from the eigenvalue spectrum, then PCA-compress.
% cc -A  : SVD over the full k-space volume (equivalent to the previous PCA approach)
% cc -M  : return the compression matrix rather than compressed data, so the
%          same matrix can be applied consistently to cal and EPI via ccapply.
ksp_gre   = bart('whiten -n', ksp_gre, ksp_noise);

% Keep the minimum number of PCA components whose cumulative explained variance
% reaches cfg.cc_energy_thresh. Lower bound: 2R for SENSE feasibility.
ksp_mat  = reshape(permute(ksp_gre, [4 1 2 3]), Ncoils, []);
C        = real((ksp_mat * ksp_mat') / size(ksp_mat, 2));
evals    = sort(real(eig(C)), 'descend');
clear ksp_mat C;

cumvar      = cumsum(evals) / sum(evals);
Nvcoils_e   = find(cumvar >= cfg.cc_energy_thresh, 1, 'first');
Nvcoils     = max(Nvcoils_e, 2 * R);
Nvcoils     = min(Nvcoils, Ncoils);
fprintf('  Energy threshold %.4f: %d components needed\n', cfg.cc_energy_thresh, Nvcoils_e);
fprintf('  Lower bound 2R = %d  →  selected Nvcoils = %d\n', 2 * R, Nvcoils);

if cfg.interactive
    figure;
    yyaxis left;
    semilogy(evals, 'b.-', 'MarkerSize', 10);
    ylabel('Eigenvalue (sample covariance)');
    yyaxis right;
    plot(cumvar * 100, 'r-');
    yline(cfg.cc_energy_thresh * 100, 'r--', sprintf('%.1f%%', cfg.cc_energy_thresh * 100));
    ylabel('Cumulative explained variance (%)');
    xline(Nvcoils, 'g-', sprintf('Nvcoils = %d', Nvcoils));
    if Nvcoils > Nvcoils_e
        xline(Nvcoils_e, 'b--', sprintf('energy → %d', Nvcoils_e));
    end
    xlabel('Component index');
    title(sprintf('Coil compression eigenvalue spectrum  (selected Nvcoils = %d)', Nvcoils));
end
clear Nvcoils_e;

cc_matrix = bart(sprintf('cc -p %d -A -M', Nvcoils), ksp_gre);
ksp_gre   = bart(sprintf('ccapply -p %d', Nvcoils), ksp_gre, cc_matrix);
save(strrep(cfg.fn.gre, '.h5', '.mat'), 'ksp_gre', '-v7.3');

%% STEP 3 — Sensitivity maps
% Done here so ksp_gre is freed BEFORE the large EPI data is loaded.
% Peak memory: ksp_noise + ksp_gre + smaps_raw  (no EPI yet).
fn_smaps = fullfile(cfg.datdir, sprintf('recon/smaps_%s.mat', cfg.SENSEmethod));

if cfg.doSENSE
    fn_smaps_valid = false;
    if exist(fn_smaps, 'file')
        tmp = load(fn_smaps, 'Nvcoils');
        fn_smaps_valid = isfield(tmp, 'Nvcoils') && tmp.Nvcoils == Nvcoils;
        if ~fn_smaps_valid
            fprintf('  Smaps file Nvcoils mismatch — recomputing.\n');
        end
        clear tmp;
    end

    if fn_smaps_valid
        fprintf('Loading precomputed sensitivity maps from %s\n', fn_smaps);
        load(fn_smaps, 'smaps');
    else
        fprintf('Estimating sensitivity maps via %s...\n', cfg.SENSEmethod);
        tic
            [smaps_raw, emaps] = makeSmaps(ksp_gre, cfg.SENSEmethod);
        toc
        smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
            Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, ...
            cfg.SENSEmethod, cfg.threshold_mask);
        save(fn_smaps, 'smaps_raw', 'emaps', 'smaps', 'Nvcoils', '-v7.3');
    end
    clear smaps_raw emaps fn_smaps_valid;

    % Uncomment if x-direction alignment between GRE and EPI is needed:
    % smaps = flip(smaps, 1);
end

if cfg.interactive
    interactive4D(abs(smaps));
end

% ksp_gre is no longer needed — free it before loading EPI
clear ksp_gre;

%% STEP 4 — Calibration data  →  odd/even phase offsets a, trajectory kxo/kxe
% Small compared to EPI; process fully before loading EPI.
fprintf('Loading calibration data...\n');
try
    ksp_cal_raw = single(orc_read(cfg.fn.cal));
catch ME
    error('preprocess: Failed to read calibration file ''%s''.\n  %s', cfg.fn.cal, ME.message);
end

assert(Nfid == size(ksp_cal_raw, 1), ...
    'preprocess: Calibration Nfid (%d) does not match noise Nfid (%d).', ...
    size(ksp_cal_raw, 1), Nfid);
fprintf('  Max |Re(cal)|: %g,  Max |Im(cal)|: %g\n', ...
    max(real(ksp_cal_raw(:))), max(imag(ksp_cal_raw(:))));

% Whiten and coil-compress calibration data.
% Do NOT squeeze after whitening: keeping the singleton in dim 3 places coils
% in MATLAB dim 4 (= BART dim 3), which is where ccapply expects them.
ksp_cal = bart('whiten -n', ...
    reshape(permute(ksp_cal_raw, [1 3 2]), Nfid, [], 1, Ncoils), ksp_noise);
clear ksp_cal_raw;

ksp_cal = squeeze(bart(sprintf('ccapply -p %d', Nvcoils), ksp_cal, cc_matrix));
ksp_cal = permute(ksp_cal, [1 3 2]);  % → [Nfid, Nvcoils, N_cal]

% Load pre-computed k-space trajectory (cycles/cm) and apply gradient delay
try
    load(cfg.fn.kxoe, 'kxo', 'kxe');
catch ME
    error('preprocess: Could not load kxoe file ''%s''.\n  %s', cfg.fn.kxoe, ME.message);
end
kxo = kxo / 100;
kxe = kxe / 100;

fprintf('Applying k-space center offset: %.2f samples\n', cfg.delay);
kxo = interp1(1:Nfid, kxo, (1:Nfid) - 0.5 - cfg.delay, 'linear', 'extrap');
kxe = interp1(1:Nfid, kxe, (1:Nfid) - 0.5 - cfg.delay, 'linear', 'extrap');

% Reshape calibration data and estimate odd/even phase offsets
% [a(1) = constant offset (rad), a(2) = linear term (rad/fov)]
ksp_cal  = reshape(permute(ksp_cal, [1 3 2]), Nfid, ETL, [], Nvcoils);
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
    error('preprocess: Could not load sampling log ''%s''.\n  %s', cfg.fn.samp_log, ME.message);
end
[Nframes, Nshots, ETL, ~] = size(schedules);

%% STEP 6 — EPI data  (frame-by-frame to bound memory)
fprintf('Opening EPI archive: %s\n', cfg.fn.epi);
try
    epi_archive = GERecon('Archive.Load', cfg.fn.epi);
catch ME
    error('preprocess: Failed to open EPI archive ''%s''.\n  %s', cfg.fn.epi, ME.message);
end

shots_per_frame       = ETL * Nshots;
total_shots_expected  = shots_per_frame * Nframes;
if epi_archive.FrameCount < total_shots_expected
    error('preprocess: EPI archive has %d shots but %d expected (%d frames × %d shots/frame).', ...
        epi_archive.FrameCount, total_shots_expected, Nframes, shots_per_frame);
end

% Pre-allocate the output file so each frame can be written immediately,
% avoiding any need to hold the full time-series in memory.
fprintf('Pre-allocating output file: %s\n', cfg.fn.recon);
mf = matfile(cfg.fn.recon, 'Writable', true);
mf.ksp_epi_zf = complex(zeros(Nx, Ny, Nz, Nvcoils, Nframes, 'single'));

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

    % ── Coil-compress → [Nfid, Nvcoils, shots_per_frame] ─────────────────
    ksp_frame = squeeze(bart(sprintf('ccapply -p %d', Nvcoils), ksp_frame, cc_matrix));
    ksp_frame = permute(ksp_frame, [1 3 2]);

    % ── Reshape → [Nfid, ETL, Nshots, Nvcoils] ────────────────────────────
    % Matches the layout the original produced via reshape+permute on the
    % full dataset before splitting across frames.
    ksp_frame = reshape(ksp_frame, Nfid, Nvcoils, ETL, Nshots);
    ksp_frame = permute(ksp_frame, [1 3 4 2]);

    % ── Grid along kx ──────────────────────────────────────────────────────
    ksp_frame_cart = hmriutils.epi.rampsampepi2cart( ...
        ksp_frame, kxo, kxe, Nx, fov(1)*100, 'nufft');

    % ── Odd/even phase correction ───────────────────────────────────────────
    ksp_frame_cart = hmriutils.epi.epiphasecorrect(ksp_frame_cart, a);

    % ── Scatter readouts into zero-filled volume ───────────────────────────
    ksp_frame_zf = zeros(Nx, Ny, Nz, Nvcoils, 'single');
    for shot_idx = 1:Nshots
        for echo = 1:ETL
            iy = schedules(frame, shot_idx, echo, 1);
            iz = schedules(frame, shot_idx, echo, 2);

            if any(ksp_frame_zf(:, iy, iz, :) ~= 0)
                warning('preprocess: Overwriting frame %d, ky=%d, kz=%d. Check schedule.', ...
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

imgs_mc = zeros(Nx, Ny, Nz, Nvcoils, NtestFrames);
for frame = 1:NtestFrames
    imgs_mc(:, :, :, :, frame) = toppe.utils.ift3(ksp_test(:, :, :, :, frame));
end

% Coil combination
if cfg.doSENSE
    img_final = squeeze(sum(imgs_mc .* conj(smaps), 4));  % Matched-filter combination
else
    img_final = squeeze(sqrt(sum(abs(imgs_mc).^2, 4)));   % Root-sum-of-squares
end

if cfg.interactive
    interactive4D(abs(permute(img_final,   [2 3 1 4])));
    interactive4D(angle(permute(img_final, [2 3 1 4])));
end
