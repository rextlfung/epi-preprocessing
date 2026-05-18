%{
run_preprocessing.m — Batch driver for the 3D EPI preprocessing pipeline.

Loads session-level configuration from config.m, then calls preprocess()
for each sequence in cfg.seqnames.  To process a single sequence, set
cfg.seqnames = {'caipi_ts'} (one entry).

Edit config.m to change sequences or parameters.
%}

run('./config.m');   % Loads session-level cfg (seqnames, datdir, gre, params...)

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    seqname = cfg.seqnames{i};
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), seqname);
    cfg_seq = set_seq_paths(cfg, seqname);
    try
        preprocess(cfg_seq);
    catch ME
        fprintf('ERROR in ''%s'': %s\nContinuing...\n', seqname, ME.message);
    end
end

fprintf('\nBatch complete.\n');
