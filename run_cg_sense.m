%{
run_cg_sense.m — Batch CG-SENSE reconstruction of all temporal frames.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame via the built-in cg_sense(). Output is saved as a .mat
file containing the complex image and all reconstruction parameters, for each
sequence in cfg.seqnames. No external toolboxes beyond standard MATLAB are
required.

Set cfg.interactive = false (the default) when running in batch — graphics
cannot be displayed from workers.

Dependencies: cg_sense.m, recon_frames.m
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), cfg.seqnames{i});
    cfg_seq  = set_seq_paths(cfg, cfg.seqnames{i});

    num_iter = cfg_seq.num_iter;  % captured by the anonymous function below
    datdir   = strcat(cfg_seq.datdir, 'recon/');
    fn_recon = fullfile(datdir, sprintf('%s_recon_cgs_i%d.mat', cfg_seq.seqname, cfg_seq.num_iter));

    try
        [img, seq_params, runtime_s] = recon_frames(cfg_seq, '', @(data, smaps) cg_sense(data, smaps, num_iter));

        Nx = seq_params.Nx;  Ny = seq_params.Ny;  Nz = seq_params.Nz;
        fov = seq_params.fov;  volumeTR = seq_params.volumeTR;
        Nvcoils = seq_params.Nvcoils;  Nframes = seq_params.Nframes;

        SENSEmethod      = cfg_seq.SENSEmethod;
        threshold_mask   = cfg_seq.threshold_mask;
        doSENSE          = cfg_seq.doSENSE;
        cc_energy_thresh = cfg_seq.cc_energy_thresh;
        seqname          = cfg_seq.seqname;

        fprintf('Saving reconstruction to %s\n', fn_recon);
        save(fn_recon, ...
            'img', 'num_iter', ...
            'SENSEmethod', 'threshold_mask', 'doSENSE', 'cc_energy_thresh', ...
            'seqname', 'Nx', 'Ny', 'Nz', 'fov', 'volumeTR', 'Nvcoils', 'Nframes', ...
            'runtime_s', '-v7.3');
    catch ME
        fprintf('ERROR [%s]: %s\nSkipping...\n', cfg_seq.seqname, ME.message);
    end
end

fprintf('\nBatch complete.\n');
