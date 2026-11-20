# Explicit-bath data and case settings

## Generic zero-temperature spin–boson cases

[generic_zero_T](generic_zero_T/) contains the 16 cases formed from four spectra and four requested mode counts \(N_b=20,40,60,80\). Frequencies are in units of \(\omega_c=1\), with \(\Delta=0.1\), zero bias, \(V=\sigma_z\), and a factorized localized spin and bare bath vacuum as the target initial state. The design and propagation window is \(T=150/\omega_c\), with 1,501 output times.

| Spectrum | Selected parameters | Hard numerical interval in units of \(\omega_c\) |
| --- | --- | --- |
| sub-Ohmic | \(s=0.25\), \(\alpha=0.03\) | \(10^{-7}\) to 4 |
| Ohmic | \(s=1\), \(\alpha=0.03\) | \(10^{-7}\) to 4 |
| super-Ohmic | \(s=2\), \(\alpha=0.03\) | \(10^{-7}\) to 4 |
| TLS-shaped | Lorentzian amplitude 0.010, center 0.1, half-width 0.015, Ohmic background strength 0.004 | \(10^{-5}\) to 4 |

The archived production driver used mD₂ multiplicity 8, seed 1, spread branches, amplitude and displacement perturbation scales \(10^{-8}\) and \(10^{-6}\), a TDVP shift of \(10^{-10}\), and `ode113` tolerances \(10^{-7}\) and \(10^{-9}\) with maximum step \(0.05/\omega_c\). Representative cache settings are recorded in [actual_cached_settings.json](generic_zero_T/solver_checks/actual_cached_settings.json). The LMD anchor uses \(n_A=3\), \(c_A=6\), and \(p=1\) for sub-Ohmic, Ohmic, and TLS, with \(p=2\) for super-Ohmic.

- `curves/` contains one trajectory per method and case. `P_HEOM` is the population comparison column for the aggregate. `P_reference`, where present, is the original fixed BSDO–mD₂ working trajectory.
- `population_summary.csv` gives non-ID population metadata and errors. `canonical_id_generic_summary.csv` supplies ID. The ID construction uses 500 time samples, 2,000 candidate frequencies, and 30,000 points for the continuum target. See [canonical ID settings](generic_zero_T/canonical_id/run_manifest.json).
- `bcf_N20_N40_N80.csv` and `bcf_N60.csv` supply the original BCF scores for the other seven methods. `bcf_modes/` stores the oscillator lists for those BCF benchmarks and `bcf_settings/` records their grids. The BCF and population batches have different output grids.
- In the archived CSVs, `MDM` is the centroid mode-density variant (MDM-c). It uses equal masses of \(J(\omega)d\omega/\omega\) and a \(J\)-weighted centroid in each interval. Population ranks, mean ranks, and Spearman correlations use the canonical ID rows.

## Bacteriochlorophyll-a

[BChl_zero_T](BChl_zero_T/) contains the zero-temperature dimer with \(H_S=J_{12}\sigma_x\), \(V=\sigma_z\), and a 200 fs propagation window. The representative case uses \(J_{12}=363\,\mathrm{cm}^{-1}\), half-width \(8\,\mathrm{cm}^{-1}\), and a 0.1 to \(4000\,\mathrm{cm}^{-1}\) spectral interval. The main fixed comparison uses a 60-mode BSDO bath and multiplicity 8. See [reference_gate.json](BChl_zero_T/solver_checks/reference_gate.json) for solver settings and convergence checks. Its multiplicity check failed. The saved trajectory serves as a fixed comparison.

- `bchl_population_curves.csv` contains `BSDO_reference` and the non-ID curves at requested mode counts 10, 20, 30, and 40. `curves/` and `modes/` contain the matching canonical ID runs.
- `spectrum_and_BCF/` contains the broadened spectrum, the continuum BCF, mode tables, and construction settings. The spectral axes use wavenumbers and the propagation uses time in fs. The constructor converts the spectrum to angular-frequency units.
- `sensitivity/` contains case-specific reference settings and curves for the broadening and electronic-coupling grid. The 18-case aggregate uses LMD-A/B/S at 20 and 40 modes with linewidths 4, 8, and \(16\,\mathrm{cm}^{-1}\) and couplings 250, 363, and \(500\,\mathrm{cm}^{-1}\). `bchl_panel_width_by_panel.csv` and its summary record the width statistics. A separate `response_reorganization/` directory gives the response and high-frequency moment records.

## Finite temperature, controls, and robustness

| Location | What to use |
| --- | --- |
| [finite_temperature](finite_temperature/) | The October comparison for Ohmic and TLS at \(\beta\omega_c=5\), \(\Delta/\omega_c=0.1\), \(T\omega_c=150\), and total signed-frequency mode budgets 40 and 80. The selected mD₂ multiplicities are 8 for Ohmic and 10 for TLS. `thermal_main_curves.csv` and `thermal_main_errors.csv` are the assembled comparison. `run_inputs/` contains the matching LMD/ID records and `MDM-c/` the physical MDM-c then thermofield-doubling records. Per-method JSONs record numerical choices. |
| [controls](controls/) | Kernel–anchor ablation, pure-dephasing checks, exponent sensitivity, finite-bath ED comparisons, timing, and mD₂ convergence tables. The original manifests give the case settings where applicable. |
| [multilevel](multilevel/) | Independent-bath and FMO channel and mode tables with their model metadata. These are construction and response records. |
| [robustness](robustness/) | The November diagnostic batch. `job_plan.csv` gives each condition, coupling scale, cutoff, initial spin state, design and propagation windows, and multiplicity. Each explicit-bath job folder contains its trajectory, modes, and `status.json` with the exact configuration. `paired_checks.csv` gives within-method differences. Reference-based errors appear in `observable_errors.csv` only for admitted observables and windows. |
