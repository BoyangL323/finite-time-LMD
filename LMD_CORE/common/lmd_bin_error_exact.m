function [Eb, wbar, lam2] = lmd_bin_error_exact(wgrid, Jgrid, i0, i1, tgrid, wfun, mode, cfg)
%LMD_BIN_ERROR_EXACT  Exact one-point local error for bin [wgrid(i0), wgrid(i1)].
%
%   [Eb, wbar, lam2] = lmd_bin_error_exact(..., mode, cfg)
%
%   M0     = int_Ib J dw ;  M1 = int_Ib w J dw ;  wbar = M1/M0 ;  lam2 = M0
%
% There are TWO error functionals, one per variational target, and they must
% not be mixed (mode selects which):
%
%   mode = 'l2'    (default, ROUTE A)   target = || Delta C ||^2_{L2[0,T]}
%       phi_b(t) = (1/M0) int_Ib J(w) exp(-i(w - wbar) t) dw
%       Eb       = M0^2 int_0^T w(t) |1 - phi_b(t)|^2 dt
%       wfun is the time weight; the route-A target uses w(t) = 1.
%
%   mode = 'phase' (ROUTE B)            target = influence phase at T
%       Eb = | int_Ib J(w) [ G_T(w) - G_T(wbar) ] dw |^2
%          = | int_0^T (T - t) Delta C_b(t) dt |^2
%       This is the modulus squared of ONE complex number, i.e. it is COHERENT
%       in time.  Feeding wfun = @(t)(T-t) into the 'l2' branch instead would
%       give int (T-t)|Delta C_b|^2 dt, a weighted BCF L2 norm, which is a
%       DIFFERENT functional and is not the route-B target.  wfun is ignored
%       here; T is taken from max(tgrid).
%
%   mode = 'phase_sys' (LMD-S CONTROLLED DIAGNOSTIC)
%       Eb = sum_q p_q |int_Ib J(w)[G_T(w-Omega_q)-
%              G_T(wbar-Omega_q)]dw|^2.
%       Omega_q and p_q are the Bohr frequencies and normalized squared
%       coupling matrix elements of the configured (H_S,V).  Its narrow-bin
%       kernel is exactly lmd_kernel_sysaware; LMD-S is not a default route.
%
% Consistency of each mode with its own seed, in the narrow-bin limit:
%   'l2'    Eb -> M0^2 (sigma^2/2)^2 int w(t) t^4 dt        ->  dmu ~ Jt^(1/3)
%   'phase' Eb -> M0^2 (sigma^2/2)^2 |G_T''(wbar)|^2        ->  dmu ~ [Jt|G_T''|]^(1/3)
%   'phase_sys' replaces |G_T''|^2 by sum_q p_q|G_T''(wbar-Omega_q)|^2.
% test_lmd_bin_error_exact evaluates convergence to these limits at runtime;
% no historical, unversioned measurement is used as an acceptance result.
%
% The 'phase' branch is evaluated as a single integral of G_T(w) - G_T(wbar)
% rather than as the difference of two O(M0 T^2) numbers, which removes the
% leading cancellation.
%
% This is the exact one-point replacement for the two historical heuristics:
%   H_b = M0 (T sigma_eff)^2  is the Taylor limit of Eb for dw_b*T << 1;
%   q_b = T sigma_b           is the same shape proxy with the M0 weight removed.
% Neither is used any more.  The exact expression saturates for wide bins
% instead of inheriting the Taylor proxy's unbounded t^4 penalty.
%
% RESET THE CACHE with  lmd_bin_error_exact('reset').  The cache is keyed on a
% signature of (wgrid, Jgrid, tgrid, wfun) and is rebuilt automatically when the
% signature changes, so a stale cache cannot leak between spectra.

persistent CACHE

if nargin == 1 && (ischar(wgrid) || isstring(wgrid)) && strcmpi(char(wgrid), 'reset')
    CACHE = [];
    Eb = []; wbar = []; lam2 = [];
    return;
end
if nargin < 6 || isempty(wfun), wfun = @(t) ones(size(t)); end
if nargin < 7 || isempty(mode), mode = 'l2'; end
if nargin < 8, cfg = []; end
mode = lower(char(mode));

wgrid = wgrid(:);  Jgrid = Jgrid(:);  tgrid = tgrid(:);
if numel(wgrid) ~= numel(Jgrid) || numel(wgrid) < 2
    error('lmd_bin_error_exact:GridSize', ...
        'wgrid and Jgrid must have equal length >= 2.');
end
if any(~isfinite(wgrid)) || any(diff(wgrid) <= 0) || any(~isfinite(Jgrid))
    error('lmd_bin_error_exact:BadGrid', ...
        'wgrid must be finite/strictly increasing and Jgrid must be finite.');
end
if numel(tgrid) < 2 || any(~isfinite(tgrid)) || any(diff(tgrid) < 0)
    error('lmd_bin_error_exact:BadTimeGrid', ...
        'tgrid must contain at least two finite nondecreasing points.');
end
i0 = min(numel(wgrid), max(1, round(i0)));
i1 = min(numel(wgrid), max(1, round(i1)));
if i1 < i0
    error('lmd_bin_error_exact:BadBin', 'Require i1 >= i0.');
end

% Do not key a scientific cache by endpoint/sum checksums: two different
% spectra can share those summaries.  Compare every input that can affect the
% result exactly.  Evaluating wfun also captures values hidden in an anonymous
% function closure that func2str alone would miss.
wprobe = wfun(tgrid);
wprobe = wprobe(:);
if numel(wprobe) == 1, wprobe = wprobe*ones(size(tgrid)); end
if numel(wprobe) ~= numel(tgrid) || any(~isfinite(wprobe))
    error('lmd_bin_error_exact:BadTimeWeight', ...
        'wfun(tgrid) must return one finite value per time or a finite scalar.');
end
cfgsig = local_cfg_signature(cfg);
same_inputs = ~isempty(CACHE) ...
    && isequaln(CACHE.wgrid,wgrid) && isequaln(CACHE.Jgrid,Jgrid) ...
    && isequaln(CACHE.tgrid,tgrid) && isequaln(CACHE.wprobe,wprobe) ...
    && strcmp(CACHE.mode,mode) && isequaln(CACHE.cfgsig,cfgsig);
if ~same_inputs
    CACHE = struct('wgrid',wgrid,'Jgrid',Jgrid,'tgrid',tgrid, ...
                   'wprobe',wprobe,'mode',mode,'cfgsig',cfgsig, ...
                   'map',containers.Map('KeyType','char','ValueType','any'));
end
key = sprintf('%d_%d', i0, i1);
if isKey(CACHE.map, key)
    v = CACHE.map(key);
    Eb = v(1); wbar = v(2); lam2 = v(3);
    return;
end

if i1 <= i0
    wbar = wgrid(max(i0, 1));
    lam2 = 0.0;
    Eb   = 0.0;
    CACHE.map(key) = [Eb, wbar, lam2];
    return;
end

ww = wgrid(i0:i1);
JJ = Jgrid(i0:i1);
dw = lmd_composite_weights(ww);

M0 = sum(JJ .* dw);
M1 = sum(ww .* JJ .* dw);
if ~(M0 > 0)
    wbar = sqrt(ww(1) * ww(end));
    lam2 = 0.0;
    Eb   = 0.0;
    CACHE.map(key) = [Eb, wbar, lam2];
    return;
end
wbar = min(max(M1 / M0, ww(1)), ww(end));
lam2 = M0;

switch mode
    case 'l2'
        t  = tgrid;
        nt = numel(t);
        phi = zeros(nt, 1);
        % Chunk over frequency so the kernel matrix never exceeds ~2e6 entries.
        chunk = max(1, floor(2e6 / max(nt, 1)));
        for a = 1:chunk:numel(ww)
            b = min(numel(ww), a + chunk - 1);
            phi = phi + exp(-1i * (t * (ww(a:b) - wbar).')) * (JJ(a:b) .* dw(a:b));
        end
        phi = phi / M0;

        Eb = (M0^2) * trapz(t, wprobe .* abs(1 - phi).^2);

    case 'phase'
        T  = max(tgrid);
        dG = lmd_kernel_G1(ww, T, cfg) - lmd_kernel_G1(wbar, T, cfg);
        Eb = abs(sum(JJ .* dw .* dG))^2;

    case 'phase_sys'
        if isempty(cfg)
            error('lmd_bin_error_exact:MissingCfg','phase_sys requires cfg.');
        end
        T = max(tgrid);
        [Omega,pw] = local_system_channels(cfg);
        if isempty(pw)
            Eb = 0;
        else
            pw = pw/sum(pw);
            Eb = 0;
            for q = 1:numel(pw)
                dG = lmd_kernel_G1(ww-Omega(q),T,cfg) ...
                    - lmd_kernel_G1(wbar-Omega(q),T,cfg);
                Eb = Eb + pw(q)*abs(sum(JJ.*dw.*dG))^2;
            end
        end

    otherwise
        error('lmd_bin_error_exact:BadMode', ...
              'mode must be ''l2'', ''phase'', or ''phase_sys'', got %s.', mode);
end

CACHE.map(key) = [Eb, wbar, lam2];
end

function sig = local_cfg_signature(cfg)
if isempty(cfg)
    sig = struct('default_kernel',true);
    return;
end
% Only fields used by G_T belong here. Store their
% exact values and shapes rather than a lossy scalar checksum.
f = {'kernel_series_cut','kernel_series_M','T_fit','epsilon','Delta', ...
     'sigmaz','sigmax','VV1'};
sig = struct();
for k = 1:numel(f)
    if isfield(cfg,f{k}), sig.(f{k}) = cfg.(f{k}); else, sig.(f{k}) = []; end
end
end

function [Omega,pw] = local_system_channels(cfg)
Hs = 0.5*cfg.epsilon*cfg.sigmaz - 0.5*cfg.Delta*cfg.sigmax;
Hs = (Hs+Hs')/2;
[U,D] = eig(Hs);
e = real(diag(D));
Vt = U'*cfg.VV1*U;
[ea,eb] = ndgrid(e,e);
Omega = ea(:)-eb(:);
pw = abs(Vt(:)).^2;
keep = isfinite(Omega) & isfinite(pw) & pw > 0;
Omega = Omega(keep); pw = pw(keep);
end
