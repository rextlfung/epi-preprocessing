%{
run_rss.m — Batch root-sum-of-squares (RSS) coil combination reconstruction.

Reads the zero-filled k-space produced by run_preprocessing.m and reconstructs
every temporal frame by applying a 3-D IFFT followed by RSS coil combination.
No sensitivity maps are used. Output is saved as a .mat file containing the
image and sequence parameters, for each sequence in cfg.seqnames. No external
toolboxes beyond standard MATLAB and TOPPE are required.

Dependencies: recon_frames.m, toppe.utils.ift3
%}

%% ── Configuration ────────────────────────────────────────────────────────────
run('./config.m');
for p = cfg.addpaths; addpath(p{1}); end

fprintf('Batch: %d sequence(s) in %s\n', numel(cfg.seqnames), cfg.datdir);

for i = 1:numel(cfg.seqnames)
    fprintf('\n[%d/%d] %s\n', i, numel(cfg.seqnames), cfg.seqnames{i});
    cfg_seq = set_seq_paths(cfg, cfg.seqnames{i});

    datdir   = strcat(cfg_seq.datdir, 'recon/');
    fn_recon = fullfile(datdir, sprintf('%s_recon_rss.mat', cfg_seq.seqname));

    try
        [img, seq_params, runtime_s] = recon_frames(cfg_seq, '', @(data, ~) sqrt(sum(abs(toppe.utils.ift3(data)).^2, 4)));

        Nx = seq_params.Nx;  Ny = seq_params.Ny;  Nz = seq_params.Nz;
        fov = seq_params.fov;  volumeTR = seq_params.volumeTR;
        Nvcoils = seq_params.Nvcoils;  Nframes = seq_params.Nframes;
        seqname = cfg_seq.seqname;

        fprintf('Saving reconstruction to %s\n', fn_recon);
        save(fn_recon, ...
            'img', ...
            'seqname', 'Nx', 'Ny', 'Nz', 'fov', 'volumeTR', 'Nvcoils', 'Nframes', ...
            'runtime_s', '-v7.3');
    catch ME
        fprintf('ERROR [%s]: %s\nSkipping...\n', cfg_seq.seqname, ME.message);
    end
end

fprintf('\nBatch complete.\n');
