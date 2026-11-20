function cfg = lmd_default_config()
%LMD_DEFAULT_CONFIG  Shared user-editable defaults for the revised LMD-TW1 package.
%
% ---------------------------------------------------------------------------
% COUPLING AND BCF CONVENTION (audited 2026-08-12, see LMD_convention_audit)
% ---------------------------------------------------------------------------
% Spectral density is defined WITHOUT a factor of pi:
%
%   J(w)      = sum_k lambda_k^2 delta(w - w_k)
%   C(t)      = int_0^inf J(w) exp(-i*w*t) dw
%   C_disc(t) = sum_k lambda_k^2 exp(-i*w_k*t)
%
% The Hamiltonian that is actually propagated by DavydovD2Rabi.m carries an
% explicit factor 1/2 in the interaction term:
%
%   H = H_S + sum_k w_k b_k' b_k + (1/2) V sum_k lambda_k (b_k' + b_k)
%
% Consequences that fix every numerical prefactor:
%   (i)  the operator coupling to the system is B = (1/2) sum_k lambda_k
%        (b_k' + b_k), whose correlation function is C(t)/4;
%   (ii) the bath reorganization energy is lambda_reorg = int J(w)/(4w) dw;
%   (iii) for the power-law family J(w) = 2 alpha wc^(1-s) w^s exp(-w/wc) the
%        effective Kondo parameter is alpha_eff = alpha/4, i.e. cfg.alpha=0.03
%        corresponds to alpha_eff = 0.0075.
% Readers using J(w) = pi sum_k lambda_k^2 delta(w - w_k) must rescale.
% test_convention_dephasing_C5.m verifies the factor 1/2 numerically.
%
% No free-system / no-dephasing baseline line is plotted anywhere.

cfg = struct();

% ---- system and propagation ----
cfg.wc = 1.0;
cfg.alpha = 0.03;
cfg.Delta = 0.1;
cfg.epsilon = 0.0;
cfg.MM = 8;
cfg.Ns = 2;
cfg.sigmaz = [1,0;0,-1];
cfg.sigmax = [0,1;1,0];
cfg.VV1 = cfg.sigmaz;

cfg.t_end = 150.0;
cfg.n_t = 1501;
cfg.RelTol = 1e-6;
cfg.AbsTol = 1e-8;
cfg.MaxStep = [];
cfg.rng_seed = 1;
cfg.init_A_noise = 1e-4;
cfg.init_F_noise = 1e-2;
cfg.init_branch_mode = 'single'; % 'single' preserves legacy runs; BChl uses 'spread'
cfg.initial_system_state = [1;0]; % normalized two-level amplitude; pure-dephasing test uses |+x>
cfg.tdvp_reg = 1e-10;            % scale-aware TDVP metric shift
cfg.tdvp_solver = 'regularized'; % 'regularized' or publication cross-check 'svd'
cfg.tdvp_svd_rcond = 1e-10;      % relative singular-value cutoff for 'svd'
cfg.compute_geom = false;
cfg.store_FF = false;

% ---- frequency / integration ----
cfg.wlow_default = 1e-7;
cfg.wlow_1f_tls = 1e-5;
cfg.whigh = 4.0;
cfg.Ngrid = 20000;

% Time resolution parameter used by the panel-error audit.
cfg.T_fit = 150.0;
cfg.id_Nt = 4000;

% ---- LMD options ----
% gauss      : DEFAULT since 2026-08-12. Anchor sector = n_A-point Gauss rule on
%              the modified measure J/w^p over [wlow, cA/T_fit], back-transformed
%              by lam2 = v.*x.^p.  Reproduces M_{-p} .. M_{2 n_A - 1 - p}
%              exactly (for the usual p=1: M_{-1} .. M_{2 n_A - 2}).
% fixed      : legacy. target centroids [1,3,5]/Tfit.
% saware     : legacy. power-law q_s(x)=x^(s/3)*(1+x^2)^(-1/6) anchor constants.
% shapeaware : legacy. generic background mass anchors.
% auto       : legacy dispatch; saware for power-law spectra, shapeaware otherwise.
% The legacy modes are retained so that the A0->A3 ablation can be run inside a
% single code version (constraint C-8).
cfg.lmd_anchor_mode = 'gauss';
cfg.lmd_xA = 6.0;
cfg.lmd_T_floor = 0.0;
cfg.lmd_pA = 0.10;
cfg.lmd_use_TW = true;
cfg.lmd_use_anchors = true;

% ---- JCTC revision: variational seed, exact bin error, quadrature anchors ----
% Production choices remain fixed by the derivation.  Fields that permit a
% non-production value (notably lmd_seed_exponent) exist only so the declared
% D3 ablation is reproducible and are labelled as such in every manifest.
%
% lmd_seed        : which variational target the seed measure is derived from.
%   'a'        route A, target = ||Delta C||_{L2[0,T]}   ->  dmu ~ Jt^(1/3) dw
%   'b'        route B, target = influence phase at T    ->  dmu ~ (Jt*|G_T''|)^(1/3)
%   'legacy'   historical manuscript, ad hoc W_T = 1/sqrt(1+(wT)^2)
%   'sysaware' LMD-S controlled diagnostic: RMS-shifted route-B kernel at
%              the Bohr frequencies of (H_S,V).  It is system-specific and
%              must never be presented as the reusable/default LMD route.
% The exponent 1/3 and the kernel are the solution of ONE variational problem.
% Never change the exponent when changing the kernel (D3 is verification only).
% There is intentionally no silent scientific default.  Article runners use
% explicit method names LMD-A/LMD-B/LMD-S/LMD-legacy.  A bare 'lmd' call must
% set cfg.lmd_seed itself, otherwise lmd_discretize_wanglu raises an error.
cfg.lmd_seed        = 'require_explicit';
cfg.lmd_nA          = 3;          % anchor Gauss points (MODE count, not moments)
cfg.lmd_cA          = 6.0;        % w_A = cA/T_fit ; |G_T''| half-max at wT~5.67
cfg.lmd_anchor_p    = 'auto';     % 'auto' | 1 | 2  (auto = convergence test)
cfg.lmd_refine      = 'equalmass';% 'equalmass' | 'greedy'
cfg.lmd_greedy_rank = 'permode';  % 'abs' | 'permode'
cfg.lmd_npts_bin    = 1;          % 1 = mainline ; 2 = all-bin control/upper bound
cfg.lmd_n2_budget   = 0;          % N2 in B0 + N2 = N_bg ; exactly accounted
cfg.lmd_seed_exponent = 1/3;      % production derivation; other values are D3 ablation only
cfg.kernel_series_cut = 10;       % wT switch point for lmd_kernel_G2
cfg.kernel_series_M   = 80;       % series terms below the switch

% ---- 1/f spectrum ----
cfg.j1f_amp = 2e-3;
cfg.j1f_w0 = 1e-3;
cfg.j1f_wcut = 1.0;

% ---- TLS / Lorentzian spectrum ----
cfg.tls_amp = 0.010;
cfg.tls_center = 0.1;
cfg.tls_gamma = 0.015;
cfg.tls_background_alpha = 0.004;

end
