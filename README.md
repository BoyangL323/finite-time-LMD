# Finite-Time Local-Moment Discretization of Bosonic Baths

This package accompanies the study by Boyang Liu and Yang Zhao. It contains MATLAB bath constructors, a multiple Davydov \(D_2\) example, numerical records from the study, and the rules used to analyze the trajectories.

| Directory | Contents |
| --- | --- |
| [LMD_CORE](LMD_CORE/README.md) | Separate LMD-A, LMD-B, and LMD-S entry points, shared numerical routines, the original group mD₂ propagator, and a sub-Ohmic example. |
| [DATA](DATA/README.md) | Explicit-bath case inputs and outputs under `mDA/` and hierarchy-depth trajectories and fitting settings under `HEOM/`. Selected supplementary figures are filed with their data. |
| [CALCULATIONS](CALCULATIONS/README.md) | How to align trajectories, calculate errors and window statistics, select references, and reproduce rankings. |

## Quick start

From the repository root in MATLAB, run

```matlab
addpath(fullfile(pwd, 'LMD_CORE', 'examples'))
[modes, trajectory] = demo_subohmic_lmd_md2();
```

The example uses 20 bath oscillators, a design and propagation time of \(30/\omega_c\), multiplicity 2, and the study's sub-Ohmic spectral parameters. It returns oscillator frequencies and weights, \(P_1(t)\), and the variational-state squared norm. The generic production benchmark used four spectra, \(T=150/\omega_c\), and mode budgets of 20, 40, 60, and 80. The [core README](LMD_CORE/README.md) explains the constructor inputs and coupling convention. The [case README](DATA/mDA/README.md) gives the production settings and points to per-case records.

The example calls the original group `DavydovD2Rabi.m` implementation. The archived production trajectories used a derived variant with a different placement of the TDVP regularizer. The example demonstrates how to use LMD with mD₂. To reproduce the reported analysis, use the archived CSVs and the [calculation rules](CALCULATIONS/README.md).
