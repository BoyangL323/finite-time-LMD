function [nodes, weights, info] = lmd_gauss_from_measure(x, wts, n)
%LMD_GAUSS_FROM_MEASURE  n-point Gauss rule for a discrete measure.
%
%   [nodes, weights] = lmd_gauss_from_measure(x, wts, n)
%
%   x    : quadrature abscissae (dense grid), column or row
%   wts  : measure weights at x, ALREADY including dw (i.e. sigma(x).*dx);
%          use lmd_composite_weights to build dx
%   n    : number of Gauss nodes
%
% Discrete Stieltjes (Lanczos) recurrence -> alpha, beta
% Jacobi matrix   J = diag(alpha) + diag(sqrt(beta(2:end)), +-1)
% nodes   = eig(J)
% weights = beta(1) * V(1,:).^2        (V = eigenvectors of J)
%
% The rule integrates every polynomial of degree <= 2n-1 exactly against the
% supplied discrete measure.  This is the shared primitive behind
%   * the BSDO / T-TEDOPA star-mapping baseline (measure = J(w) dw), and
%   * the low-frequency anchor sector (measure = J(w) w^{-p} dw), where the
%     back-transform lam2 = weights .* nodes.^p turns exactness in
%     {1, w, ..., w^{2n-1}} against J/w^p into exactness of the moments
%     M_{-p}, M_{1-p}, ..., M_{2n-1-p} of J itself.  The commonly used p=1
%     case is M_{-1},...,M_{2n-2}; p=2 instead ends at M_{2n-3}.
%
% Power moments + Hankel matrices are deliberately NOT used: that route is
% hopelessly ill-conditioned beyond n ~ 20.
%
% info fields:
%   alpha, beta  recurrence coefficients (beta(1) = total mass)
%   n_used       nodes actually produced (< n if the recurrence broke down)
%   beta_min     smallest beta(2:end); a tiny/negative value flags loss of
%                orthogonality, which is the expected failure mode for
%                sub-Ohmic measures at large n in double precision
%   truncated    true if n_used < n

x   = double(x(:));
wts = double(wts(:));
if numel(x) ~= numel(wts)
    error('lmd_gauss_from_measure:Size', 'x and wts must have equal length.');
end
ok  = isfinite(x) & isfinite(wts) & (wts > 0);
x   = x(ok);
wts = wts(ok);
if isempty(x)
    error('lmd_gauss_from_measure:Empty', 'Measure is empty after filtering.');
end
n = max(1, min(round(n), numel(x)));

N     = numel(x);
alpha = zeros(n, 1);
beta  = zeros(n, 1);
p0    = zeros(N, 1);
p1    = ones(N, 1);
beta(1) = sum(wts);
n_used  = n;

for k = 1:n
    nk = sum(wts .* p1 .* p1);
    if ~(nk > 0) || ~isfinite(nk)
        n_used = k - 1;
        break;
    end
    alpha(k) = sum(wts .* x .* p1 .* p1) / nk;
    if k + 1 <= n
        if k > 1
            p2 = (x - alpha(k)) .* p1 - beta(k) * p0;
        else
            p2 = (x - alpha(k)) .* p1;
        end
        bk = sum(wts .* p2 .* p2) / nk;
        if ~(bk > 0) || ~isfinite(bk)
            n_used = k;                 % recurrence exhausted the measure
            break;
        end
        beta(k+1) = bk;
        p0 = p1;
        p1 = p2;
    end
end

n_used = max(1, min(n_used, n));
a = alpha(1:n_used);
b = beta(1:n_used);

if n_used == 1
    nodes   = a(1);
    weights = b(1);
else
    off = sqrt(b(2:end));
    Jm  = diag(a) + diag(off, 1) + diag(off, -1);
    [V, D] = eig((Jm + Jm.') / 2);
    nodes   = real(diag(D));
    weights = b(1) * (V(1, :).').^2;
    [nodes, idx] = sort(nodes, 'ascend');
    weights = weights(idx);
end

nodes   = nodes(:);
weights = max(weights(:), 0.0);

info = struct();
info.alpha     = a;
info.beta      = b;
info.n_used    = n_used;
info.truncated = (n_used < n);
if numel(b) > 1
    info.beta_min = min(b(2:end));
else
    info.beta_min = NaN;
end
info.mass      = beta(1);
info.mass_rule = sum(weights);
end
