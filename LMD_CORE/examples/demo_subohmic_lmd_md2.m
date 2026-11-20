function [modes, trajectory] = demo_subohmic_lmd_md2()
%DEMO_SUBOHMIC_LMD_MD2 Construct a sub-Ohmic LMD-S bath and propagate it.

coreRoot = fileparts(fileparts(mfilename('fullpath')));
originalPath = path;
restorePath = onCleanup(@() path(originalPath));
addpath(fullfile(coreRoot, 'common'), fullfile(coreRoot, 'LMD-S'), ...
    fullfile(coreRoot, 'mDA'));

cfg = lmd_default_config();
cfg.MM = 2;
cfg.T_fit = 30;
cfg.t_end = cfg.T_fit;
cfg.n_t = 301;
cfg.MaxStep = 0.1;
modeBudget = 20;
spec = lmd_get_spectrum('subohmic', cfg);
[omega, lambda, info] = lmd_s(spec, modeBudget, cfg);
modes = table(omega, lambda.^2, info.is_anchor(:), ...
    'VariableNames', {'omega', 'lambda_squared', 'is_anchor'});
trajectory = propagate_md2(omega, lambda, cfg);
fprintf('Sub-Ohmic LMD-S: %d modes, %d samples, final P1 = %.6g\n', ...
    height(modes), height(trajectory), trajectory.P1(end));
end
