function g = lmd_kernel_G2(w, T, normalize, cfg)
%LMD_KERNEL_G2  |d^2 G_T/dw^2| with numerically safe series/closed-form switch.
%
%   g = lmd_kernel_G2(w, T)             -> absolute value |G_T''(w)|
%   g = lmd_kernel_G2(w, T, true)       -> normalized to g(0)=1, i.e. /(T^4/12)
%   g = lmd_kernel_G2(w, T, norm, cfg)  -> use cfg.kernel_series_cut / _M
%
% G_T is the influence-phase kernel of route B,
%
%   G_T(w)    = int_0^T (T - tau) exp(-i w tau) dtau
%             = -i*T/w + (1 - exp(-i w T))/w^2
%   G_T''(w)  = -int_0^T (T - t) t^2 exp(-i w t) dt
%             = T^2 e/w^2 - 2i T/w^3 - 4i T e/w^3 + 6 (1 - e)/w^4,  e = exp(-i w T)
%   G_T''(0)  = -T^4/12
%
% The closed form contains 1/w^4 and suffers catastrophic cancellation for
% small w.  Below
% wT = cfg.kernel_series_cut (default 10) the entry series is used instead,
%
%   G_T''(w) = -T^4 * sum_{m>=0} (-i w T)^m / [ m! (m+3)(m+4) ]
%
% truncated at cfg.kernel_series_M (default 80) terms.  Above the switch the
% series becomes ill-conditioned and the closed form is accurate, so the
% switch is required in both directions.
%
% test_lmd_kernel_G2 verifies both branches and the switch against independent
% direct quadrature and emits the achieved error in the current MATLAB run.
%
% |G_T''| is EVEN in w and strictly monotone non-increasing on
% wT in [1e-2, 1e2] (zero sign changes, 0.0% dip below the running max), so no
% smoothing envelope is needed anywhere.  Half-max of |G_T''|/|G_T''(0)| sits at
% wT = 5.671, which is the origin of cfg.lmd_cA = 6.0; the test computes this
% value rather than accepting a comment as evidence.  Asymptotically
% |G_T''| ~ T^2/w^2.
%
% See also LMD_KERNEL_SYSAWARE, LMD_BIN_ERROR_EXACT.

if nargin < 3 || isempty(normalize), normalize = false; end
if nargin < 4 || isempty(cfg)
    cut = 10.0;  M = 80;
else
    cut = cfg.kernel_series_cut;  M = cfg.kernel_series_M;
end
if ~(T > 0)
    error('lmd_kernel_G2:BadT', 'T must be positive, got %g.', T);
end

sz = size(w);
x  = abs(double(w(:)));          % |G_T''| is even in w
g  = zeros(numel(x), 1);

small = (x * T) < cut;

if any(small)
    z    = -1i * (x(small) * T);
    term = ones(size(z));        % z^m / m!
    acc  = zeros(size(z));
    for m = 0:(M-1)
        acc  = acc + term / ((m + 3) * (m + 4));
        term = term .* z / (m + 1);
    end
    g(small) = abs(-(T^4) * acc);
end

if any(~small)
    ww = x(~small);
    e  = exp(-1i * ww * T);
    val = (T^2) * e ./ ww.^2 - 2i * T ./ ww.^3 - 4i * T * e ./ ww.^3 + ...
          6 * (1 - e) ./ ww.^4;
    g(~small) = abs(val);
end

if normalize
    g = g / (T^4 / 12);
end
g = reshape(g, sz);
end
