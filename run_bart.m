%{
run_bart.m — Batch L1-SENSE reconstruction via BART pics.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame via BART pics. Output is saved as a .mat file containing
the complex image and all reconstruction parameters, for each sequence in
cfg.seqnames.

Dependencies: BART (https://mrirecon.github.io/bart/), recon_frames.m
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), cfg.seqnames{i});
    cfg_seq = set_seq_paths(cfg, cfg.seqnames{i});

    bart_cmd = sprintf('pics -R W:7:0:%g -R T:7:0:%g -i %d -S', cfg_seq.lamb_l1, cfg_seq.lamb_tv, cfg_seq.num_iter);
    fprintf('BART command: %s\n', bart_cmd);

    datdir   = strcat(cfg_seq.datdir, 'recon/');
    fn_recon = fullfile(datdir, sprintf('%s_recon_bart_l1_r%.4f_tv_r%.4f.mat', cfg_seq.seqname, cfg_seq.lamb_l1, cfg_seq.lamb_tv));

    try
        [img, seq_params, runtime_s] = recon_frames(cfg_seq, '', @(data, smaps) bart(bart_cmd, data, smaps));

        Nx = seq_params.Nx;  Ny = seq_params.Ny;  Nz = seq_params.Nz;
        fov = seq_params.fov;  volumeTR = seq_params.volumeTR;
        Nvcoils = seq_params.Nvcoils;  Nframes = seq_params.Nframes;

        lamb_l1          = cfg_seq.lamb_l1;
        lamb_tv          = cfg_seq.lamb_tv;
        num_iter         = cfg_seq.num_iter;
        SENSEmethod      = cfg_seq.SENSEmethod;
        threshold_mask   = cfg_seq.threshold_mask;
        doSENSE          = cfg_seq.doSENSE;
        cc_energy_thresh = cfg_seq.cc_energy_thresh;
        seqname          = cfg_seq.seqname;

        fprintf('Saving reconstruction to %s\n', fn_recon);
        save(fn_recon, ...
            'img', 'bart_cmd', ...
            'lamb_l1', 'lamb_tv', 'num_iter', ...
            'SENSEmethod', 'threshold_mask', 'doSENSE', 'cc_energy_thresh', ...
            'seqname', 'Nx', 'Ny', 'Nz', 'fov', 'volumeTR', 'Nvcoils', 'Nframes', ...
            'runtime_s', '-v7.3');
    catch ME
        fprintf('ERROR [%s]: %s\nSkipping...\n', cfg_seq.seqname, ME.message);
    end
end

fprintf('\nBatch complete.\n');
