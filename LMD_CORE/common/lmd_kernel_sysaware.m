function [K, info] = lmd_kernel_sysaware(w, cfg)
%LMD_KERNEL_SYSAWARE  System-aware route-B kernel (LMD-S diagnostic only).
%
%   K(w) = sqrt(sum_q p_q |G_T''(w-Omega_q)|^2),
%   p_q  = |<a|V|b>|^2 / sum_ab |<a|V|b>|^2,
%   Omega_q = E_a-E_b.
%
% K is normalized to unit maximum on the supplied grid.  The same channel
% weights and shifts are used by lmd_bin_error_exact(...,'phase_sys',cfg), so
% the seed measure and refinement functional are one controlled construction.
% This route depends on H_S and V and is therefore never the reusable default.

if nargin < 2 || isempty(cfg), cfg = lmd_default_config(); end

Hs = 0.5*cfg.epsilon*cfg.sigmaz - 0.5*cfg.Delta*cfg.sigmax;
Hs = (Hs+Hs')/2;
[U,D] = eig(Hs);
e = real(diag(D));
Vt = U'*cfg.VV1*U;
P = abs(Vt).^2;
[ea,eb] = ndgrid(e,e);
Omega = ea(:)-eb(:);
pw = P(:);
keep = isfinite(Omega) & isfinite(pw) & pw > 0;
Omega = Omega(keep);
pw = pw(keep);

sz = size(w);
x = double(w(:));
if isempty(pw)
    K = reshape(ones(size(x)),sz);
    info = struct('Omega',[],'p',[],'eigenvalues',e,'p_zero',0, ...
        'note','V has no nonzero energy-basis element; kernel set to one');
    return;
end
pw = pw/sum(pw);
K = zeros(size(x));
for q = 1:numel(pw)
    gq = lmd_kernel_G2(x-Omega(q),cfg.T_fit,true,cfg);
    K = K + pw(q)*gq.^2;
end
K = sqrt(max(K,0));
mx = max(K);
if mx > 0, K = K/mx; else, K(:) = 1; end
K = reshape(K,sz);

info = struct('Omega',Omega,'p',pw,'eigenvalues',e, ...
    'p_zero',sum(pw(abs(Omega)<1e-12)), ...
    'note','RMS shifted-channel kernel; normalized only by a global constant');
end
