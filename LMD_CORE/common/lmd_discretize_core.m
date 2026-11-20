function [omega, lambda, meta] = lmd_discretize_core(spec, Nb, cfg, route)
%LMD_DISCRETIZE_CORE Shared LMD branch extracted from the archived solver.
% Call one of the public LMD-A, LMD-B or LMD-S entry points.

if isempty(cfg), cfg = lmd_default_config(); end
defcfg = lmd_default_config();
def_fields = fieldnames(defcfg);
for ii_def = 1:numel(def_fields)
    if ~isfield(cfg, def_fields{ii_def})
        cfg.(def_fields{ii_def}) = defcfg.(def_fields{ii_def});
    end
end
wlow = spec.wlow;
whigh = spec.whigh;
wgrid = exp(linspace(log(wlow), log(whigh), cfg.Ngrid)).';
Jgrid = spec.Jfun(wgrid);
Jgrid = max(Jgrid(:), 0.0);
C0 = trapz(wgrid, Jgrid);
if ~(C0 > 0)
    error('Nonpositive C0 for spectrum %s.', spec.name);
end

[omega, lam2, meta] = local_lmd_tw1(spec, Nb, cfg, wgrid, Jgrid, C0);

omega = omega(:); lam2 = lam2(:);
keep = isfinite(omega) & isfinite(lam2) & omega > 0 & lam2 > 0;
omega = omega(keep); lam2 = lam2(keep);
[omega, idx] = sort(omega, 'ascend');
lam2 = lam2(idx);
lambda = sqrt(max(lam2, 0.0));
% keep the per-mode annotations aligned with the surviving, sorted modes
for fn = {'is_anchor','rule'}
    f = fn{1};
    if isfield(meta, f) && numel(meta.(f)) == numel(keep)
        v = meta.(f); v = v(keep); meta.(f) = v(idx);
    end
end

meta.method = lower(route);
meta.label = route;
meta.Nb_requested = Nb;
meta.Nb_target = Nb;                       % C-6: always report both
meta.Nb_eff = numel(omega);                % C-6: nonzero propagated modes
meta.C0 = C0;
meta.wlow = wlow;
meta.whigh = whigh;
meta.wgrid = wgrid;
meta.Jgrid = Jgrid;
meta.spectrum = spec.name;
meta.cfg_snapshot = local_cfg_snapshot(cfg);
meta.coverage = local_coverage_stats(omega, lam2, cfg, whigh);   % C-7
end

function s = local_cfg_snapshot(cfg)
% Everything that changes the discretization, for the run manifest (section 5).
s = struct();
f = {'T_fit','whigh','Ngrid','lmd_seed','lmd_anchor_mode','lmd_nA','lmd_cA', ...
     'lmd_anchor_p','lmd_refine','lmd_greedy_rank','lmd_npts_bin','lmd_n2_budget', ...
     'lmd_seed_exponent', ...
     'lmd_use_TW','lmd_use_anchors','lmd_xA','lmd_pA','kernel_series_cut', ...
     'kernel_series_M','alpha','wc','Delta','epsilon'};
for k = 1:numel(f)
    if isfield(cfg, f{k}), s.(f{k}) = cfg.(f{k}); end
end
end

function c = local_coverage_stats(omega, lam2, cfg, whigh)
% C-7: representative-frequency coverage and the per-band mode census.
wA = cfg.lmd_cA / cfg.T_fit;
if whigh > 2
    % Historical dimensionless five-spectrum bands (kept byte-for-byte in
    % meaning for the generic regression suite).
    bands = [0, wA; wA, 1; 1, 2; 2, whigh];
else
    % Molecular inputs are commonly in atomic units, where whigh << 1.  The
    % historical [1,2] edges would then create reversed/empty bands.  Retain the
    % physical anchor boundary and add monotone fractional high-frequency bands.
    be = unique(sort([0, min(max(wA,0),whigh), 0.25*whigh, ...
                      0.50*whigh, 0.75*whigh, whigh]));
    be = be([true, diff(be) > max(eps(whigh),1e-15)]);
    bands = [be(1:end-1).', be(2:end).'];
end
c = struct();
c.omega_max = max([omega(:); -Inf]);
c.omega_min = min([omega(:);  Inf]);
c.omega_max_over_whigh = c.omega_max / whigh;
c.band_edges  = bands;
c.band_counts = zeros(size(bands,1), 1);
c.band_mass   = zeros(size(bands,1), 1);
for b = 1:size(bands,1)
    m = omega >= bands(b,1) & omega < bands(b,2);
    if b == size(bands,1), m = omega >= bands(b,1) & omega <= bands(b,2); end
    c.band_counts(b) = sum(m);
    c.band_mass(b)   = sum(lam2(m));
end
if numel(omega) > 1
    d = diff(sort(omega));
    d = d(d > 0);
    if isempty(d)
        c.min_spacing = 0; c.t_recurrence = Inf;
    else
        c.min_spacing  = min(d);
        c.t_recurrence = 2*pi / min(d);      % t_rec = 2 pi / min_k dw_k
    end
else
    c.min_spacing = NaN; c.t_recurrence = Inf;
end
c.t_recurrence_over_Tfit = c.t_recurrence / cfg.T_fit;
c.max_lambda_over_omega  = max(sqrt(max(lam2,0)) ./ max(omega, eps));
c.any_negative_weight    = any(lam2 < 0);
c.any_out_of_band        = any(omega < 0) || any(omega > whigh * (1 + 1e-9));
end

function [omega, lam2, meta] = local_lmd_tw1(spec, Nb, cfg, wgrid, Jgrid, C0)
%LOCAL_LMD_TW1  LMD core: seed measure -> anchor sector -> background rules.
%
% Structure (spec sections 3.3 - 3.5):
%   1. seed density  mu = [ Jt * K ]^(1/3)   with K fixed by cfg.lmd_seed
%   2. anchor sector [wlow, w_A]: n_A-point Gauss rule on J/w^p (mode 'gauss'),
%      or one of the legacy anchor heuristics
%   3. background sector [w_A, whigh]: equal-mass or E_b-greedy bins, each
%      carrying a 1-point (centroid) rule, optionally upgraded to 2 points with
%      exact budget accounting  B0 + N2 = N_bg
%
% The exponent 1/3 and the kernel K are the solution of ONE variational problem
% (equal error per bin). Never retune the exponent when the kernel changes.

wgrid = wgrid(:);
Jgrid = Jgrid(:);
Jt    = Jgrid / max(C0, 1e-300);
Tfit  = cfg.T_fit;

% ---- 1. seed measure --------------------------------------------------------
[mu, seed_meta] = local_lmd_seed(wgrid, Jt, cfg);
Wenv = seed_meta.kernel;

% ---- 2. anchor sector -------------------------------------------------------
anchor_mode = lower(char(cfg.lmd_anchor_mode));
if strcmp(anchor_mode, 'auto')
    if strcmp(spec.kind, 'powerlaw')
        anchor_mode = 'saware';
    else
        anchor_mode = 'shapeaware';
    end
end

omega_a = [];
lam2t_a = [];
anchor_meta = struct('mode', anchor_mode, 'p', NaN, 'wA', NaN, 'nA', 0, 'info', []);

if cfg.lmd_use_anchors && Nb >= 4 && strcmp(anchor_mode, 'gauss')
    % ---- anchor sector: local Gauss rule on the modified measure J/w^p ----
    wA = cfg.lmd_cA / Tfit;
    wA = min(max(wA, wgrid(1) * 10), wgrid(end) * 0.5);
    nA = max(1, min(cfg.lmd_nA, Nb - 1));
    [p, pinfo] = lmd_select_anchor_power(spec, cfg);
    % The anchor measure must terminate exactly at w_A, not at the last grid
    % point below it: otherwise the sliver [wgrid(k), w_A] is counted by neither
    % sector and M0 acquires a spurious O(J(w_A) dw) defect.
    m   = (wgrid >= wgrid(1)) & (wgrid <= wA);
    wa  = wgrid(m);
    Ja  = Jt(m);
    if wa(end) < wA
        wa(end+1,1) = wA;
        Ja(end+1,1) = interp1(wgrid, Jt, wA, 'linear', 'extrap');
    end
    sig = Ja ./ wa.^p;
    dwa = lmd_composite_weights(wa);
    [x, v, ginfo] = lmd_gauss_from_measure(wa, sig .* dwa, nA);
    omega_a = x;
    lam2t_a = v .* x.^p;                     % <-- back-transform, do not omit
    left    = wA;
    anchor_meta.p     = p;
    anchor_meta.wA    = wA;
    anchor_meta.nA    = numel(omega_a);
    anchor_meta.info  = ginfo;
    anchor_meta.power_selection = pinfo;
    edges_anchor = [wgrid(1); wA];
else
    % ---- legacy anchor heuristics (kept for the A0->A3 ablation, C-8) ----
    [edges_anchor, anchors_legacy, left] = local_legacy_anchor_edges(spec, Nb, cfg, wgrid, Jt, mu, anchor_mode);
    anchor_meta.legacy_targets = anchors_legacy;
    if numel(edges_anchor) > 1
        [omega_a, lam2t_a] = local_bin_moments_grid(wgrid, Jt, edges_anchor);
        anchor_meta.nA = numel(omega_a);
    end
    anchor_meta.wA = left;
end

% ---- 3. background sector ---------------------------------------------------
nA_used = numel(omega_a);
N_bg    = max(Nb - nA_used, 1);

% Exact 2-point budget accounting (constraint C-1).  B0 bins are laid down, N2
% of them are upgraded to a 2-point rule, so the mode count is exactly
%       B0 * 1 + N2 * 1  =  B0 + N2  =  N_bg .
% Nothing is ever "free": every extra node is debited from the bin count.
if cfg.lmd_npts_bin >= 2
    B0 = max(1, ceil(N_bg / 2));             % all-bin 2pt: CONTROL / upper bound
    N2 = N_bg - B0;                          % odd N_bg leaves one 1-point bin
else
    N2 = max(0, min(round(cfg.lmd_n2_budget), floor(N_bg / 2)));
    B0 = N_bg - N2;                          % selective upgrade
end

i_lo = find(wgrid >= left, 1, 'first');
if isempty(i_lo), i_lo = 1; end
i_lo = max(1, min(i_lo, numel(wgrid) - 1));
i_hi = numel(wgrid);

tgrid = linspace(0, Tfit, max(400, min(1600, cfg.id_Nt))).';
wfun  = seed_meta.time_weight;
emode = seed_meta.err_mode;

greedy_info = [];
switch lower(char(cfg.lmd_refine))
    case 'equalmass'
        edges_bg = local_equal_mass_edges(wgrid, mu, left, wgrid(i_hi), B0);
    case 'greedy'
        [edges_bg, greedy_info] = lmd_greedy_refine(wgrid, Jgrid, B0, tgrid, wfun, ...
            min(B0, 4), mu, cfg.lmd_greedy_rank, i_lo, i_hi, emode, cfg);
    otherwise
        error('lmd_discretize_wanglu:BadRefine', ...
              'cfg.lmd_refine must be ''equalmass'' or ''greedy'', got %s.', cfg.lmd_refine);
end
% the background sector must start exactly where the anchor sector ended
edges_bg = local_monotone_edges(edges_bg(:), left, wgrid(i_hi));

% per-bin point counts: 1 everywhere, then upgrade the N2 worst bins to 2
nbg  = numel(edges_bg) - 1;
npts = ones(nbg, 1);
if N2 > 0
    Ebin = zeros(nbg, 1);
    for b = 1:nbg
        jb0 = local_index_of(wgrid, edges_bg(b));
        jb1 = local_index_of(wgrid, edges_bg(b+1));
        Ebin(b) = lmd_bin_error_exact(wgrid, Jgrid, jb0, jb1, tgrid, wfun, emode, cfg);
    end
    [~, ord] = sort(Ebin, 'descend');
    npts(ord(1:min(N2, nbg))) = 2;
end

[omega_b, lam2t_b, rule_b] = local_bin_rules_grid(wgrid, Jt, edges_bg, npts);

% ---- assemble ---------------------------------------------------------------
omega = [omega_a(:); omega_b(:)];
lam2  = C0 * [lam2t_a(:); lam2t_b(:)];
is_anchor = [true(numel(omega_a),1); false(numel(omega_b),1)];

meta = struct();
meta.variant       = sprintf(['LMD: seed=%s, anchor=%s, refine=%s, ' ...
                              'background 1pt/2pt = %d/%d'], ...
                              seed_meta.name, anchor_mode, lower(char(cfg.lmd_refine)), ...
                              sum(npts == 1), sum(npts == 2));
meta.seed          = seed_meta;
meta.anchor        = anchor_meta;
meta.anchor_mode   = anchor_mode;
meta.anchors       = omega_a(:);
meta.edges_anchor  = edges_anchor(:);
meta.edges         = [edges_anchor(1:end-1); edges_bg(:)];
meta.edges_bg      = edges_bg(:);
meta.npts_bin      = npts(:);
meta.rule          = [repmat({'gauss-anchor'}, numel(omega_a), 1); rule_b(:)];
meta.is_anchor     = is_anchor;
meta.mu            = mu;
meta.Wenv          = Wenv;
meta.Jtilde        = Jt;
meta.Tfit          = Tfit;
meta.greedy        = greedy_info;
meta.budget        = struct('Nb_target', Nb, 'n_anchor', nA_used, 'N_bg', N_bg, ...
                            'B0', B0, 'N2', N2, 'n_modes_built', numel(omega), ...
                            'accounting_ok', (B0 + N2) == N_bg);
end

% ---------------------------------------------------------------------------
function [mu, info] = local_lmd_seed(wgrid, Jt, cfg)
%LOCAL_LMD_SEED  Seed density of the three (four) variational targets.
%
%   'a'        route A: target = ||Delta C_b||^2_{L2[0,T]}. Because |e^{-i wbar t}|
%              = 1 and int_0^T t^4 dt = T^5/5 carries no w dependence, equal error
%              per bin gives  dmu ~ Jt^(1/3) dw  with NO time-window envelope.
%   'b'        route B: target = unmodulated endpoint influence-phase surrogate.
%              Per-bin contribution
%              ~ C0 M0 (sigma^2/2) |G_T''(wbar)|, so  dmu ~ [Jt |G_T''|]^(1/3) dw.
%              The envelope APPEARS, it is not appended.
%   'legacy'   historical manuscript, ad hoc W_T = 1/sqrt(1 + (wT)^2).
%   'sysaware' LMD-S: route B shifted to the Bohr-frequency channels of H_S,V.
%
% cfg.lmd_use_TW = false switches the kernel off entirely (ablation D1); it has
% no effect on route 'a', which has no kernel by construction.
% The SEED and the per-bin ERROR FUNCTIONAL must come from the same target.
% info.err_mode / info.time_weight are what lmd_bin_error_exact is called with:
%   route A / legacy -> mode 'l2'       with w(t) = 1
%   route B          -> mode 'phase'     (one coherent endpoint defect)
%   LMD-S            -> mode 'phase_sys' (weighted channel defects)
% Using mode 'l2' with the weight (T-t) for route B would silently optimize
% int (T-t)|Delta C|^2 dt, a weighted BCF norm, not the influence phase.
Tfit = cfg.T_fit;
seed = lower(char(cfg.lmd_seed));
switch seed
    case 'a'
        K = ones(size(wgrid));
        wfun = @(t) ones(size(t));
        emode = 'l2';
        nm = 'route A  [Jt]^(1/3), no envelope';
    case 'b'
        K = lmd_kernel_G2(wgrid, Tfit, true, cfg);
        wfun = @(t) ones(size(t));
        emode = 'phase';
        nm = 'route B  [Jt |G_T'''']]^(1/3)';
    case 'legacy'
        K = 1 ./ sqrt(1 + (wgrid * Tfit).^2);
        wfun = @(t) ones(size(t));
        emode = 'l2';
        nm = 'legacy   [Jt W_T]^(1/3)';
    case {'sysaware','s'}
        [K, channel_info] = lmd_kernel_sysaware(wgrid, cfg);
        wfun = @(t) ones(size(t));
        emode = 'phase_sys';
        nm = 'LMD-S    [Jt K_sys]^(1/3), controlled diagnostic';
    otherwise
        error('lmd_discretize_wanglu:BadSeed', ...
              ['A bare LMD call has no scientific default. Set cfg.lmd_seed ' ...
               'to ''a'', ''b'', ''sysaware'', or ''legacy''; got %s.'], seed);
end
if ~cfg.lmd_use_TW
    K = ones(size(wgrid));
    % A kernel-disabled phase seed cannot retain a phase refinement: that
    % would mix Route-A nodes with a Route-B error functional.  This ablation
    % is therefore the complete no-kernel L2 control.
    emode = 'l2';
    wfun = @(t) ones(size(t));
    nm = [nm ' (kernel disabled: full no-kernel L2 control)'];
end
pexp=cfg.lmd_seed_exponent;
if ~(isscalar(pexp)&&isfinite(pexp)&&pexp>0)
    error('lmd_discretize_wanglu:BadSeedExponent','lmd_seed_exponent must be positive.');
end
mu = max(Jt .* K, 1e-300).^pexp;
info = struct('name', nm, 'key', seed, 'kernel', K, 'exponent', pexp, ...
              'time_weight', wfun, 'err_mode', emode);
if exist('channel_info','var'), info.system_channels = channel_info; end
end

% ---------------------------------------------------------------------------
function [edges, anchors, left] = local_legacy_anchor_edges(spec, Nb, cfg, wgrid, Jt, mu, anchor_mode)
edges = wgrid(1);
left = wgrid(1);
anchors = [];
Tfit = cfg.T_fit;
if ~(cfg.lmd_use_anchors && Nb >= 4)
    return;
end
    switch anchor_mode
        case 'fixed'
            raw_targets = [1,3,5] / Tfit;
            for k = 1:numel(raw_targets)
                [edges, anchors, left] = local_append_anchor_by_centroid(edges, anchors, left, wgrid(end), raw_targets(k), wgrid, Jt, Nb);
            end
        case 'saware'
            if ~strcmp(spec.kind, 'powerlaw') || isnan(spec.s)
                raw_targets = [1,3,5] / Tfit;
            else
                Tref = max(Tfit, cfg.lmd_T_floor);
                if ~(Tref > 0), Tref = Tfit; end
                c = local_saware_anchor_constants(spec.s, cfg.lmd_xA, 3);
                raw_targets = c(:).' / Tref;
            end
            for k = 1:numel(raw_targets)
                [edges, anchors, left] = local_append_anchor_by_centroid(edges, anchors, left, wgrid(end), raw_targets(k), wgrid, Jt, Nb);
            end
        case 'shapeaware'
            [raw_targets, eguess] = local_shapeaware_anchor_targets(wgrid, mu, cfg.lmd_pA, 3);
            for k = 1:numel(raw_targets)
                if numel(anchors) >= Nb - 1, break; end
                edge_k = eguess(k+1);
                if edge_k > left*(1+1e-8) && edge_k < wgrid(end)*(1-1e-10)
                    edges(end+1,1) = edge_k; %#ok<AGROW>
                    anchors(end+1,1) = raw_targets(k); %#ok<AGROW>
                    left = edge_k;
                end
            end
        case 'none'
            % no anchors
        otherwise
            error('Unknown lmd_anchor_mode: %s', cfg.lmd_anchor_mode);
    end
edges = edges(:);
if edges(end) < left
    edges(end+1,1) = left;
end
end

function [edges, anchors, left] = local_append_anchor_by_centroid(edges, anchors, left, whigh, target, wgrid, Jt, Nb)
if numel(anchors) >= Nb - 1, return; end
if target <= left*(1+1e-8) || target >= 0.98*whigh, return; end
edge_k = local_find_centroid_edge(wgrid, Jt, left, whigh, target);
if edge_k > left*(1+1e-8) && edge_k < whigh*(1-1e-10)
    edges(end+1,1) = edge_k; %#ok<AGROW>
    anchors(end+1,1) = target; %#ok<AGROW>
    left = edge_k;
end
end

function c = local_saware_anchor_constants(s, xA, n)
x = linspace(0, xA, 20000).';
x(1) = max(1e-12, xA*1e-12);
q = (x.^(s/3)) .* ((1 + x.^2).^(-1/6));
Q = cumtrapz(x, q);
tot = Q(end);
xb = interp1(Q, x, linspace(0, tot, n+1).', 'linear', 'extrap');
xb(1) = x(1); xb(end) = xA;
c = zeros(n,1);
for i = 1:n
    mask = x >= xb(i) & x <= xb(i+1);
    xx = x(mask); qq = q(mask);
    if numel(xx) < 2
        c(i) = 0.5*(xb(i)+xb(i+1));
    else
        c(i) = trapz(xx, xx.*qq) / max(trapz(xx, qq), 1e-300);
    end
end
end

function [targets, edges] = local_shapeaware_anchor_targets(wgrid, mu, pA, n)
w = wgrid(:); q = mu(:);
Q = cumtrapz(w, q);
tot = Q(end);
QA = max(min(pA, 0.90), 1e-4) * tot;
wA = interp1(Q, w, QA, 'linear', 'extrap');
mask = w <= wA;
ww = w(mask); qq = q(mask);
if numel(ww) < n+2
    ww = w(1:min(numel(w), n+2)); qq = q(1:numel(ww));
end
QQ = cumtrapz(ww, qq); totA = QQ(end);
edges = interp1(QQ, ww, linspace(0, totA, n+1).', 'linear', 'extrap');
edges(1) = w(1); edges(end) = wA;
targets = zeros(n,1);
for i = 1:n
    maski = ww >= edges(i) & ww <= edges(i+1);
    xx = ww(maski); qqi = qq(maski);
    targets(i) = trapz(xx, xx.*qqi) / max(trapz(xx, qqi), 1e-300);
end
end

function edges = local_equal_mass_edges(x, mu, left, right, B)
if B < 1
    edges = [left; right]; return;
end
mask = x >= left & x <= right;
xx = x(mask); mm = mu(mask);
if numel(xx) < 2
    edges = linspace(left, right, B+1).'; return;
end
if xx(1) > left
    xx = [left; xx]; mm = [interp1(x, mu, left, 'linear', 'extrap'); mm];
end
if xx(end) < right
    xx = [xx; right]; mm = [mm; interp1(x, mu, right, 'linear', 'extrap')];
end
cum = cumtrapz(xx, mm); tot = cum(end);
if ~(tot > 0)
    edges = linspace(left, right, B+1).'; return;
end
keep = [true; diff(cum) > 1e-300];
cu = cum(keep); xu = xx(keep);
edges = interp1(cu, xu, linspace(0, tot, B+1).', 'linear', 'extrap');
edges(1) = left; edges(end) = right;
end

function edges = local_monotone_edges(edges, wlow, whigh)
edges = edges(:); edges(1) = wlow; edges(end) = whigh;
for i = 2:numel(edges)
    if edges(i) <= edges(i-1)
        edges(i) = edges(i-1) * (1 + 1e-12);
    end
end
edges(end) = whigh;
end

function wR = local_find_centroid_edge(w, Jt, wL, wMax, target)
mask = w >= wL & w <= wMax;
ww = w(mask); JJ = Jt(mask);
if ww(1) > wL
    ww = [wL; ww]; JJ = [interp1(w, Jt, wL, 'linear', 'extrap'); JJ];
end
cum0 = cumtrapz(ww, JJ);
cum1 = cumtrapz(ww, ww .* JJ);
cent = cum1 ./ max(cum0, 1e-300);
idx = find(cent >= target, 1, 'first');
if isempty(idx)
    wR = min(wMax*(1-1e-8), max(wL*(1+1e-8), 1.2*target));
elseif idx == 1
    wR = ww(1);
else
    wR = interp1(cent(idx-1:idx), ww(idx-1:idx), target, 'linear', 'extrap');
end
wR = min(max(wR, wL*(1+1e-8)), wMax*(1-1e-8));
end

function [omega, lam2, rule] = local_bin_rules_grid(wgrid, Jgrid, edges, npts)
%LOCAL_BIN_RULES_GRID  Per-bin local quadrature rules with mixed point counts.
%
% npts(b) = 1  -> centroid rule: omega = M1/M0, lam2 = M0.  Exact for M0 and M1
%                 of that bin, hence exact for the GLOBAL M0 and M1 of any
%                 partition.  Local order eps ~ N_bin^-2.
% npts(b) = n  -> n-point Gauss rule w.r.t. J restricted to the bin.  Exact for
%                 the first 2n bin moments; local order eps ~ N_bin^-2n.
%
% n = 1 is the mainline (constraint C-1).  n = 2 on every bin is a control /
% upper bound and must be labelled as such; selective upgrades must debit the
% budget exactly (B0 + N2 = N_bg), which the caller enforces.
B = numel(edges) - 1;
if nargin < 4 || isempty(npts), npts = ones(B,1); end
npts = npts(:);
omega = []; lam2 = []; rule = {};
for b = 1:B
    n = max(1, round(npts(min(b, numel(npts)))));
    if n <= 1
        [m0, m1] = local_moments_interval(wgrid, Jgrid, edges(b), edges(b+1));
        if m0 > 0
            w1 = min(max(m1/m0, edges(b)), edges(b+1));
        else
            w1 = sqrt(edges(b)*edges(b+1));
        end
        omega(end+1,1) = w1;             %#ok<AGROW>
        lam2(end+1,1)  = max(m0, 0.0);   %#ok<AGROW>
        rule{end+1,1}  = '1pt-centroid'; %#ok<AGROW>
    else
        [ww, JJ] = local_bin_samples(wgrid, Jgrid, edges(b), edges(b+1), 4*n);
        dw = lmd_composite_weights(ww);
        if sum(JJ .* dw) <= 0
            omega(end+1,1) = sqrt(edges(b)*edges(b+1)); %#ok<AGROW>
            lam2(end+1,1)  = 0.0;                        %#ok<AGROW>
            rule{end+1,1}  = '1pt-empty';                %#ok<AGROW>
            continue;
        end
        [x, v] = lmd_gauss_from_measure(ww, JJ .* dw, n);
        omega = [omega; x(:)];            %#ok<AGROW>
        lam2  = [lam2;  v(:)];            %#ok<AGROW>
        for q = 1:numel(x)
            rule{end+1,1} = sprintf('%dpt-gauss', n); %#ok<AGROW>
        end
    end
end
omega = omega(:); lam2 = lam2(:); rule = rule(:);
end

function [ww, JJ] = local_bin_samples(wgrid, Jgrid, wL, wR, nmin)
%LOCAL_BIN_SAMPLES  Samples of J inside [wL, wR], refined to at least nmin points.
mask = wgrid >= wL & wgrid <= wR;
ww = wgrid(mask); JJ = Jgrid(mask);
if isempty(ww) || ww(1) > wL
    ww = [wL; ww]; JJ = [interp1(wgrid, Jgrid, wL, 'linear', 'extrap'); JJ];
end
if ww(end) < wR
    ww = [ww; wR]; JJ = [JJ; interp1(wgrid, Jgrid, wR, 'linear', 'extrap')];
end
if numel(ww) < nmin
    wf = linspace(wL, wR, nmin).';
    JJ = interp1(ww, JJ, wf, 'linear', 'extrap');
    ww = wf;
end
ww = ww(:); JJ = max(JJ(:), 0.0);
end

function i = local_index_of(wgrid, w)
[~, i] = min(abs(wgrid(:) - w));
end

function [omega, lam2] = local_bin_moments_grid(wgrid, Jgrid, edges)
B = numel(edges) - 1;
omega = zeros(B,1); lam2 = zeros(B,1);
for b = 1:B
    [m0, m1] = local_moments_interval(wgrid, Jgrid, edges(b), edges(b+1));
    lam2(b) = max(m0, 0.0);
    if m0 > 0
        omega(b) = min(max(m1/m0, edges(b)), edges(b+1));
    else
        omega(b) = sqrt(edges(b)*edges(b+1));
    end
end
end

function [m0, m1] = local_moments_interval(w, J, wL, wR)
mask = w >= wL & w <= wR;
ww = w(mask); JJ = J(mask);
if isempty(ww)
    wm = sqrt(wL*wR);
    Jm = interp1(w, J, wm, 'linear', 'extrap');
    dw = wR - wL;
    m0 = Jm * dw; m1 = Jm * wm * dw; return;
end
if ww(1) > wL
    ww = [wL; ww]; JJ = [interp1(w, J, wL, 'linear', 'extrap'); JJ];
end
if ww(end) < wR
    ww = [ww; wR]; JJ = [JJ; interp1(w, J, wR, 'linear', 'extrap')];
end
m0 = trapz(ww, JJ);
m1 = trapz(ww, ww .* JJ);
end
