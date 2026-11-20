function dw = lmd_composite_weights(x)
%LMD_COMPOSITE_WEIGHTS  Trapezoidal dw weights of a (possibly non-uniform) grid.
%
%   dw = lmd_composite_weights(x)  satisfies  sum(dw .* f) == trapz(x, f)
%   for any f sampled on x.  Used to turn a sampled density sigma(x) into a
%   discrete measure sigma(x).*dw suitable for lmd_gauss_from_measure.

x = x(:);
n = numel(x);
dw = zeros(n, 1);
if n == 0
    return;
elseif n == 1
    return;                       % a single point carries no measure
end
dw(1)   = 0.5 * (x(2) - x(1));
dw(end) = 0.5 * (x(end) - x(end-1));
if n > 2
    dw(2:end-1) = 0.5 * (x(3:end) - x(1:end-2));
end
end
