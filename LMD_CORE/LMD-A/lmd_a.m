function [omega, lambda, meta] = lmd_a(spec, Nb, cfg)
%LMD_A Bath-only background kernel.
if nargin < 3 || isempty(cfg)
    cfg = lmd_default_config();
end
cfg.lmd_seed = 'a';
[omega, lambda, meta] = lmd_discretize_core(spec, Nb, cfg, 'LMD-A');
end
