function [omega, lambda, meta] = lmd_s(spec, Nb, cfg)
%LMD_S System transition-weighted finite-time response kernel.
if nargin < 3 || isempty(cfg)
    cfg = lmd_default_config();
end
cfg.lmd_seed = 'sysaware';
[omega, lambda, meta] = lmd_discretize_core(spec, Nb, cfg, 'LMD-S');
end
