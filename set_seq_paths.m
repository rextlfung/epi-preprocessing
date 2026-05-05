function cfg = set_seq_paths(cfg, seqname)
%SET_SEQ_PATHS  Populate cfg.seqname and all per-sequence cfg.fn.* fields.
%
%   cfg = set_seq_paths(cfg, seqname)
%
%   Called once per sequence in the batch loop of run_preprocessing.m, or directly when
%   running preprocess() on a single sequence.  cfg.fn.kxoe is initialised
%   with a placeholder (kxoe90.mat); preprocess() updates it after loading
%   params.m and Nx is known.
%
%   Inputs
%   ------
%   cfg      Session-level config struct (from config.m)
%   seqname  Single sequence name string, e.g. 'caipi_ts'
%
%   Output
%   ------
%   cfg      Same struct with cfg.seqname and per-sequence cfg.fn.* populated

cfg.seqname     = seqname;
cfg.seqdir      = fullfile(cfg.datdir, sprintf('seqs/%s/', seqname));
cfg.fn.params   = fullfile(cfg.datdir, sprintf('seqs/%s/params.m',           seqname));
cfg.fn.cal      = fullfile(cfg.datdir, sprintf('scanarchives/%s_cal.h5',     seqname));
cfg.fn.noise    = fullfile(cfg.datdir, sprintf('scanarchives/%s_noise.h5',   seqname));
cfg.fn.epi      = fullfile(cfg.datdir, sprintf('scanarchives/%s_epi.h5',     seqname));
cfg.fn.kxoe     = fullfile(cfg.datdir, sprintf('seqs/%s/kxoe%d.mat',         seqname, 90));
cfg.fn.samp_log = fullfile(cfg.datdir, sprintf('seqs/%s/samp_locs.mat',      seqname));
cfg.fn.recon    = fullfile(cfg.datdir, sprintf('recon/%s_epi_zf.mat',        seqname));
