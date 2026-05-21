%{
run_rss.m — Batch root-sum-of-squares (RSS) coil combination reconstruction.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame by applying a 3-D IFFT followed by RSS coil combination.
No sensitivity maps are used. Output is written to a NIfTI file for each
sequence in cfg.seqnames. No external toolboxes beyond standard MATLAB and
TOPPE are required.

Dependencies: recon_frames.m, toppe.utils.ift3
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), cfg.seqnames{i});
    cfg_seq = set_seq_paths(cfg, cfg.seqnames{i});

    datdir   = strcat(cfg_seq.datdir, 'recon/');
    fn_recon = fullfile(datdir, sprintf('%s_recon_rss.nii', cfg_seq.seqname));

    recon_frames(cfg_seq, fn_recon, @(data, ~) sqrt(sum(abs(toppe.utils.ift3(data)).^2, 4)));
end

fprintf('\nBatch complete.\n');
