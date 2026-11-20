function g = lmd_kernel_G1(w, T, cfg)
%LMD_KERNEL_G1  Complex influence-phase kernel G_T(w), safe at small w.
%
%   g = lmd_kernel_G1(w, T)        -> G_T(w)  (complex, NOT an absolute value)
%   g = lmd_kernel_G1(w, T, cfg)   -> use cfg.kernel_series_cut / _M
%
%   G_T(w) = int_0^T (T - tau) exp(-i w tau) dtau
%          = -i*T/w + (1 - exp(-i w T)) / w^2
%   G_T(0) = T^2/2
%
% This is the kernel of the route-B target: the UNMODULATED endpoint double
% integral (constant-filter / pure-dephasing influence-phase surrogate) is
% dPhi(T) = int dw drho(w) G_T(w), so the route-B per-bin error is the modulus
% squared of ONE complex number,
%
%   E_b^(Phi) = | int_{I_b} J(w) [ G_T(w) - G_T(wbar_b) ] dw |^2 ,
%
% It is not the full path-dependent influence functional of a tunnelling system,
% and it is NOT a time integral of |Delta C_b(t)|^2 with weight (T - t). The two differ:
% the first is coherent in time, the second is not.  See lmd_bin_error_exact
% with mode = 'phase'.
%
% The closed form contains 1/w^2 and cancels catastrophically for small w, so below
% cfg.kernel_series_cut (default 10) the entry series is used:
%
%   G_T(w) = T^2 * sum_{m>=0} (-i w T)^m / [ m! (m+1)(m+2) ]
%
% truncated at cfg.kernel_series_M (default 80).  Above the switch the series
% becomes ill-conditioned while the closed form is accurate.  The switch and
% both branches are checked against direct quadrature by test_lmd_kernel_G2.
%
% Consistency check that motivates the whole construction: in the narrow-bin
% limit E_b^(Phi) -> M0^2 (sigma^2/2)^2 |G_T''(wbar)|^2, so equal error per bin
% gives dmu ~ [ Jt |G_T''| ]^(1/3) dw, which is exactly the route-B seed.
%
% See also LMD_KERNEL_G2, LMD_BIN_ERROR_EXACT.

if nargin < 3 || isempty(cfg)
    cut = 10.0;  M = 80;
else
    cut = cfg.kernel_series_cut;  M = cfg.kernel_series_M;
end
if ~(T > 0)
    error('lmd_kernel_G1:BadT', 'T must be positive, got %g.', T);
end

sz = size(w);
x  = double(w(:));
g  = zeros(numel(x), 1);

small = (abs(x) * T) < cut;

if any(small)
    z    = -1i * (x(small) * T);
    term = ones(size(z));
    acc  = zeros(size(z));
    for m = 0:(M-1)
        acc  = acc + term / ((m + 1) * (m + 2));
        term = term .* z / (m + 1);
    end
    g(small) = (T^2) * acc;
end

if any(~small)
    ww = x(~small);
    g(~small) = -1i * T ./ ww + (1 - exp(-1i * ww * T)) ./ ww.^2;
end

g = reshape(g, sz);
end
