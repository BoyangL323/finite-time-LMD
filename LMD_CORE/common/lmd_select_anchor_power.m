function [p, info] = lmd_select_anchor_power(spec, cfg)
%LMD_SELECT_ANCHOR_POWER  Choose the anchor modified-measure power p.
%
%   [p, info] = lmd_select_anchor_power(spec, cfg)
%
% The anchor sector runs an n_A-point Gauss rule against J(w)/w^p on
% [wlow, w_A].  p = 1 is universally safe: it pins the moments
% M_{-1}, M_0, ..., M_{2 n_A - 2}.  p = 2 would additionally pin M_{-2}, but
% M_{-2} = int J/w^2 dw DIVERGES whenever s <= 1 (logarithmically for Ohmic,
% as a power for sub-Ohmic) and for every infrared-heavy generic spectrum, so
% the rule then collapses onto the artificial wlow cutoff.
%
% 'auto' therefore runs an explicit cutoff-sensitivity test: recompute
% int_{wlow/100}^{w_A} J/w^2 dw and require that it grow by less than 5%.  The
% 5% figure is an internal numerical criterion of the convergence test, not a
% user-facing knob (constraint C-2).
%
% info fields expose the runtime ratio and integrals so every selection is
% reproducible without relying on historical numbers in source comments.

if nargin < 2 || isempty(cfg), cfg = lmd_default_config(); end

info = struct('mode', cfg.lmd_anchor_p, 'ratio', NaN, 'M2_ref', NaN, ...
              'M2_low', NaN, 'divergent', NaN, 'reason', '');

GROWTH_TOL = 1.05;      % internal criterion of the convergence test
DECADES    = 100.0;     % cutoff is lowered by this factor for the test

wA = cfg.lmd_cA / cfg.T_fit;
wA = min(max(wA, spec.wlow * 10), spec.whigh * 0.5);

% ---- explicit user override -------------------------------------------------
if isnumeric(cfg.lmd_anchor_p) && ~isempty(cfg.lmd_anchor_p)
    p = round(cfg.lmd_anchor_p);
    if ~any(p == [1 2])
        error('lmd_select_anchor_power:BadP', 'cfg.lmd_anchor_p must be 1, 2 or ''auto''.');
    end
    info.reason = 'explicit cfg.lmd_anchor_p';
    if p == 2
        [div, ratio, I1, I2] = local_divergence_test(spec, wA, DECADES, GROWTH_TOL);
        info.ratio = ratio; info.M2_ref = I1; info.M2_low = I2; info.divergent = div;
        if div
            warning('lmd_select_anchor_power:ForcedDivergentP2', ...
                ['p = 2 forced for spectrum %s, but int J/w^2 is divergent ' ...
                 '(cutoff sensitivity %.2fx). The lowest anchor node will sit ' ...
                 'on the artificial wlow boundary.'], spec.name, ratio);
        end
    end
    return;
end

% ---- automatic selection ----------------------------------------------------
if strcmpi(char(cfg.lmd_anchor_p), 'auto')
    [div, ratio, I1, I2] = local_divergence_test(spec, wA, DECADES, GROWTH_TOL);
    info.ratio = ratio; info.M2_ref = I1; info.M2_low = I2; info.divergent = div;
    if div
        p = 1;
        info.reason = sprintf('int J/w^2 divergent (%.3fx under a %gx lower cutoff)', ratio, DECADES);
    else
        p = 2;
        info.reason = sprintf('int J/w^2 convergent (%.3fx under a %gx lower cutoff)', ratio, DECADES);
    end
    return;
end

error('lmd_select_anchor_power:BadMode', ...
      'cfg.lmd_anchor_p must be 1, 2 or ''auto''.');
end

function [divergent, ratio, I1, I2] = local_divergence_test(spec, wA, decades, tol)
% Cutoff-sensitivity test for int_0 J(w)/w^2 dw.
if isfield(spec,'ir_exponent') && isfinite(spec.ir_exponent)
    divergent = ~(spec.ir_exponent > 1.0);
    if divergent, ratio = Inf; else, ratio = 1.0; end
    I1 = NaN; I2 = NaN;
    return;
end
if strcmp(spec.kind, 'powerlaw') && isfield(spec, 's') && isfinite(spec.s)
    % Analytic shortcut: int_0 w^(s-2) dw converges iff s > 1.
    divergent = ~(spec.s > 1.0);
    if divergent
        ratio = Inf;
    else
        ratio = 1.0;
    end
    I1 = NaN; I2 = NaN;
    return;
end
n  = 20000;
g1 = exp(linspace(log(spec.wlow), log(wA), n)).';
g2 = exp(linspace(log(spec.wlow / decades), log(wA), n)).';
I1 = trapz(g1, max(spec.Jfun(g1), 0) ./ g1.^2);
I2 = trapz(g2, max(spec.Jfun(g2), 0) ./ g2.^2);
if ~(I1 > 0)
    divergent = true;  ratio = Inf;  return;
end
ratio = I2 / I1;
divergent = ~(ratio < tol);
end
