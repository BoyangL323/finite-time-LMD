# MATLAB bath constructors and mD₂ example

## Entry points

| Folder | Function | Background allocation |
| --- | --- | --- |
| `LMD-A/` | `lmd_a(spec, Nb, cfg)` | \(K_A=1\) |
| `LMD-B/` | `lmd_b(spec, Nb, cfg)` | \(K_B=|G_T''|\) |
| `LMD-S/` | `lmd_s(spec, Nb, cfg)` | Root-mean-square of the shifted magnitudes \(|G_T''(\omega-\Omega_{ab})|\), weighted by the squared coupling-operator matrix elements. The constructor included here uses the two-level model. |

All three return `[omega, lambda, meta]`. The `common/` folder holds the LMD branch extracted from the study's archived `lmd_discretize_wanglu.m` and the original kernel and quadrature routines. The routes share an infrared Gaussian anchor and allocate the remaining oscillators to equal-mass intervals of \([J(\omega)K_X(\omega)]^{1/3}\). The default background rule assigns each interval its \(J\)-weighted centroid and integrated \(J\) weight.

For example, from the repository root

```matlab
root = fullfile(pwd, 'LMD_CORE');
addpath(fullfile(root, 'common'), fullfile(root, 'LMD-B'))
cfg = lmd_default_config();
cfg.T_fit = 150;
spec = lmd_get_spectrum('subohmic', cfg);
[omega, lambda, meta] = lmd_b(spec, 40, cfg);
```

| Parameter | Meaning |
| --- | --- |
| `spec` | Use `lmd_get_spectrum` for the built-in cases. A custom spectrum supplies `name`, `kind`, `wlow`, `whigh`, and a vectorized `Jfun` over positive frequencies. For analytic power-law anchor selection, also supply the exponent `s` and set `kind='powerlaw'`. |
| `Nb` | Requested total oscillator count, including the anchor. Check `meta.Nb_eff` for the number of retained positive-weight oscillators. |
| `cfg.T_fit` | Observation time used in construction. Use the reciprocal of the frequency unit. It enters the shared anchor even for LMD-A. |
| `cfg.Delta`, `cfg.epsilon`, `cfg.VV1` | Two-level tunneling, bias, and coupling matrix. The model is \(H_S=(\epsilon\sigma_z-\Delta\sigma_x)/2\). The default `VV1` is \(\sigma_z\). |
| `cfg.Ngrid` | Size of the logarithmic frequency grid used to build the weights and centroids. |
| `cfg.lmd_nA`, `cfg.lmd_cA`, `cfg.lmd_anchor_p` | Anchor node count, nominal boundary coefficient \(\omega_A=c_A/T\), and power \(p\) for the quadrature measure \(J(\omega)\omega^{-p}d\omega\). The grid bounds constrain the numerical boundary. Defaults are \(n_A=3\), \(c_A=6\), and automatic \(p\). |
| `meta.is_anchor`, `meta.edges_bg` | Anchor flags in sorted mode order and background panel boundaries. |

The convention is \(J(\omega)=\sum_k\lambda_k^2\delta(\omega-\omega_k)\). The coupling Hamiltonian is \(V\sum_k\lambda_k(b_k^\dagger+b_k)/2\). Pass the amplitudes `lambda` directly to the propagator. The per-mode quadrature weights are `lambda.^2`.

## Multiple Davydov D₂

`mDA/DavydovD2Rabi.m` is the original group implementation distributed in `for LBY`, supplied by YYY and copied byte for byte. `propagate_md2.m` passes the oscillator arrays to its nine-argument interface and collects norm-corrected observables using `observabled2rwa.m`. The original routine adds \(10^{-10}I\) to the coherent-state overlap. For the archived production data, the derived routine instead applied a scale-aware shift to the TDVP linear system and kept the physical overlap unshifted. Production settings are documented in [DATA/mDA](../DATA/mDA/README.md).

The example uses the positive-frequency sub-Ohmic spectrum \(J=2\alpha\omega_c^{1-s}\omega^s e^{-\omega/\omega_c}\), with \(s=0.25\), \(\alpha=0.03\), \(\omega_c=1\), \(\Delta=0.1\), and \(\epsilon=0\). It starts with a localized spin and vacuum-centered variational branches, with small seeded perturbations. `cfg.MM` sets the multiplicity. `cfg.t_end`, `cfg.n_t`, `RelTol`, `AbsTol`, and `MaxStep` control propagation and output. See the [root quick start](../README.md#quick-start).
