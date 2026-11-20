function [omega, lambda, meta] = lmd_b(spec, Nb, cfg)
%LMD_B Finite-time response kernel.
if nargin < 3 || isempty(cfg)
    cfg = lmd_default_config();
end
cfg.lmd_seed = 'b';
[omega, lambda, meta] = lmd_discretize_core(spec, Nb, cfg, 'LMD-B');
end
