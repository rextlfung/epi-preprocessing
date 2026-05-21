%{
run_cg_sense.m — Batch CG-SENSE reconstruction of all temporal frames.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame via the built-in cg_sense(). Output is written to a NIfTI
file for each sequence in cfg.seqnames. No external toolboxes beyond standard
MATLAB are required.

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
    fn_recon = fullfile(datdir, sprintf('%s_recon_cgs_i%d.nii', cfg_seq.seqname, cfg_seq.num_iter));

    try
        recon_frames(cfg_seq, fn_recon, @(data, smaps) cg_sense(data, smaps, num_iter));
    catch ME
        fprintf('ERROR [%s]: %s\nSkipping...\n', cfg_seq.seqname, ME.message);
    end
end

fprintf('\nBatch complete.\n');
