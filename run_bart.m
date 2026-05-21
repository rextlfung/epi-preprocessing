%{
run_bart.m — Batch L1-SENSE reconstruction via BART pics.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame via BART pics. Output is written to a NIfTI file for each
sequence in cfg.seqnames.

Dependencies: BART (https://mrirecon.github.io/bart/), recon_frames.m
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), cfg.seqnames{i});
    cfg_seq = set_seq_paths(cfg, cfg.seqnames{i});

    bart_cmd = sprintf('pics -R W:7:0:%f -R T:7:0:%f -i 100 -S', cfg_seq.lamb_l1, cfg_seq.lamb_tv);
    fprintf('BART command: %s\n', bart_cmd);

    datdir   = strcat(cfg_seq.datdir, 'recon/');
    fn_recon = fullfile(datdir, sprintf('%s_recon_bart_l1_r%.4f_tv_r%.4f.nii', cfg_seq.seqname, cfg_seq.lamb_l1, cfg_seq.lamb_tv));

    recon_frames(cfg_seq, fn_recon, @(data, smaps) bart(bart_cmd, data, smaps));
end

fprintf('\nBatch complete.\n');
