# Trajectory selection and convergence records

## Generic zero-temperature HEOM

1. Select the continuum bath with the stated spectral cutoff, system parameters, and observation window. Match the independent fit/audit frequency grids and the propagation output time grid.
2. Require two successive ESPIRA-I ranks with stable decay rates and passing BCF thresholds. The archived test limits are \(2\times10^{-3}\) for the full relative \(L_2\) and late RMS divided by full reference RMS, and \(10^{-2}\) for maximum absolute error scaled by \(|C(0)|\). Use the saved rank-selection rule for the spectrum. The confirmed member chosen can differ by case.
3. Require at least three consecutive hierarchy depths. The last adjacent pair must satisfy \(2\times10^{-3}\) in full and late population relative \(L_2\) and maximum absolute population difference, using the preceding depth for normalization. The preceding pair must already pass, or the last pair must improve in all three metrics. The final population must remain within \([-10^{-6},1+10^{-6}]\) throughout the window.
4. Use the selected reference depths and curves in [DATA/HEOM](../DATA/HEOM/README.md). The generic gate and per-depth files give the original selection decisions. Exclude regularized \(1/f\) from the HEOM-based aggregate.

An adjacent-depth change is a convergence diagnostic. Interpret close errors with the reference and mD₂ propagation checks in view.

## Explicit-bath checks

Match spectral support, mode budget, coupling convention, Hamiltonian, initial state, and recorded propagation parameters. Use the actual nonzero mode count from the metadata when reporting a run. Validity checks in the source workflow cover finite values, physical populations and norm drift. The archived reference and multiplicity gates report their own thresholds. A higher-multiplicity comparison audits production trajectories separately from the HEOM bath-fit and depth checks.

## Bacteriochlorophyll-a

Use the prescribed BSDO–mD₂ fixed trajectory as the distance target for the representative dimer. The compact HEOM fit failed its prerequisite, and the fixed mD₂ comparison failed its multiplicity check. In the sensitivity grid, choose the case-specific fixed trajectory for each linewidth and coupling. Interpret these errors as distances within the archived comparison. The mode counts and spectral moment data can be interpreted independently of these trajectories.

## Thermal and November records

For the October thermal comparison use the saved Ohmic depth-4 and selected TLS depth-5 curves. The TLS selection record and older automatic depth-4 audit are both available. The published-basis curves use the stated multiplicities and the revised MDM-c followed by thermofield doubling.

For the November robustness tests, use the reference selection for each observable and time window. Read the case's `reference_windows.csv`, join it with the job status, and preserve `depth_unresolved` and absent references. The original-condition population records reuse the designated generic reference. A run without an admitted reference can contribute within-method paired diagnostics. The ranking screen is limited to complete cases with finite depth margins.
