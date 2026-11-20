function spec = lmd_get_spectrum(name, cfg)
%LMD_GET_SPECTRUM  Five Wang-Lu zero-temperature benchmark spectra.
%
% Supported names:
%   subohmic, ohmic, superohmic, one_over_f, tls
%
% Returned spec fields:
%   name, display, kind, s, wlow, whigh, Jfun

if nargin < 2 || isempty(cfg), cfg = lmd_default_config(); end
name = lower(strrep(strtrim(name), '-', '_'));
name = strrep(name, '1/f', 'one_over_f');
name = strrep(name, '1f', 'one_over_f');

spec = struct();
spec.name = name;
spec.whigh = cfg.whigh;
spec.alpha = cfg.alpha;
spec.wc = cfg.wc;

switch name
    case 'subohmic'
        s = 0.25;
        spec.display = 'sub-Ohmic';
        spec.kind = 'powerlaw';
        spec.s = s;
        spec.wlow = cfg.wlow_default;
        spec.Jfun = @(w) J_powerlaw_exp(w, cfg.alpha, s, cfg.wc);

    case 'ohmic'
        s = 1.0;
        spec.display = 'Ohmic';
        spec.kind = 'powerlaw';
        spec.s = s;
        spec.wlow = cfg.wlow_default;
        spec.Jfun = @(w) J_powerlaw_exp(w, cfg.alpha, s, cfg.wc);

    case 'superohmic'
        s = 2.0;
        spec.display = 'super-Ohmic';
        spec.kind = 'powerlaw';
        spec.s = s;
        spec.wlow = cfg.wlow_default;
        spec.Jfun = @(w) J_powerlaw_exp(w, cfg.alpha, s, cfg.wc);

    case 'one_over_f'
        spec.display = 'regularized 1/f';
        spec.kind = 'generic';
        spec.s = NaN;
        spec.wlow = cfg.wlow_1f_tls;
        A = cfg.j1f_amp; w0 = cfg.j1f_w0; wcut = cfg.j1f_wcut;
        spec.Jfun = @(w) A ./ (w + w0) .* exp(-w ./ wcut);

    case 'tls'
        spec.display = 'TLS Lorentzian';
        spec.kind = 'generic';
        spec.s = NaN;
        spec.wlow = cfg.wlow_1f_tls;
        A = cfg.tls_amp; nu = cfg.tls_center; gam = cfg.tls_gamma;
        abg = cfg.tls_background_alpha; wc = cfg.wc;
        spec.Jfun = @(w) A * (gam.^2 ./ ((w - nu).^2 + gam.^2)) + ...
            2 * abg * w .* exp(-w ./ wc);

    otherwise
        error('Unknown spectrum name: %s', name);
end
end

function J = J_powerlaw_exp(w, alpha, s, wc)
w = max(w, 0.0);
J = 2 * alpha * (wc^(1-s)) .* (w.^s) .* exp(-w ./ wc);
end
