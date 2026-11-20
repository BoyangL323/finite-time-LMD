# HEOM trajectories and bath-fit settings

| Location | Case and reference choice |
| --- | --- |
| [generic_zero_T](generic_zero_T/) | Four finite-window references used to rank the generic population calculations. The chosen depths are 5 for sub-Ohmic and 3 for Ohmic, super-Ohmic, and TLS. The regularized \(1/f\) depth ladder is retained as a diagnostic. Read `generic_heom_gate.json` for the final qualification decision and use the `*_heom_reference.csv` curves for the four admitted cases. |
| [finite_temperature](finite_temperature/) | \(\beta\omega_c=5\) Ohmic and TLS. Use `ohmic_beta5_2571f74d9d6c/D4.csv` and `tls_beta5_depth5_b05abddb9777/D5.csv` for the selected comparison. The TLS depth-5 decision is documented in `accepted_reference.json`. The preceding depth-4 sequence and its original audit are retained. |
| [BChl_diagnostics](BChl_diagnostics/) | Zero-temperature dimer at \(J_{12}=363\,\mathrm{cm}^{-1}\) and half-width \(8\,\mathrm{cm}^{-1}\). The compact bath representation failed its prerequisite audit. The saved positive-frequency depth-1 to depth-3 trajectories are diagnostics. See `canonical_depth_gate.json` and `bchl_compact_bath_gate.json`. |
| [robustness](robustness/) | November attempts at enhanced Ohmic coupling and the baseline TLS case. Per-case input, fit and depth audits accompany `depth_*.csv`. Both attempts are recorded as `depth_unresolved` in the job status. |

## Software and fitting

The generic and BChl archive used Python with the QuTiP 5.2 series. The upstream environment declared `qutip>=5.2,<5.3`, `numpy>=1.26`, and `scipy>=1.11`. A per-run patch version is available where stated in the saved metadata. The October finite-temperature and November audit files record QuTiP **5.2.3**.

The generic comparison fits the hard-truncated continuum \(C(t)\) on its finite time grid with QuTiP's ESPIRA-I decomposition. The bath-fit audits report the fit mode, rank-selection rule, numerical fit and verification grids, exponent count, decay rates, and errors. Typically the fit and independent continuum audit use 20,000 and 40,000 frequency points, respectively. The per-depth JSON files record solver settings. Check the final gate for acceptance, since run-level JSONs describe the trajectory at the time it was generated.

The thermal correlation uses signed-frequency weights at the prescribed inverse temperature. The archived runs use the raw correlation with \(Q=\sigma_z/2\). The selected ESPIRA-I ranks are 16 and 18 per component for Ohmic and TLS, with 32 and 36 exponentials. The solver uses BDF with relative and absolute tolerances of \(10^{-8}\) and \(10^{-10}\). See the corresponding `thermal_fit.input.json`, `thermal_fit.json`, and depth inputs for details.

The BChl diagnostic folder contains two distinct records. The positive-frequency trajectory files give the historical depth comparison for its stated environment. The compact fit files document the separate finite-window shared-pole bath audit and its result. See the [trajectory rules](../../CALCULATIONS/trajectory_selection.md) for the comparison policy.
