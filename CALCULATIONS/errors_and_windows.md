# Error metrics and window statistics

## Pairing trajectory samples

For each method, pair the logged test values \(x_i\) with the comparison values \(y_i\) on the same ordered time grid. Use the case and budget identifiers in the data index. Check that time points match before computing norms. The benchmark used uniform output grids. A full or late window selects samples on both endpoints, including the midpoint in the late window.

The spin population is \(P_1=(1+\langle\sigma_z\rangle)/2\). In generic mD₂ curve CSVs, `P_method` and `P_HEOM` already contain \(P_1\). The historical generic HEOM column named `Pz` also contains \(P_1\) and is used directly. The thermal HEOM exports distinguish `P1` and the actual spin expectation `Pz`. In the BChl CSVs, pair the logged site population columns.

## Time window statistics

Define \(e_i=x_i-y_i\) and let \(I_W\) select the output samples in window \(W\).

\[
\operatorname{RMS}_W(e)=
\left(\frac{1}{|I_W|}\sum_{i\in I_W}|e_i|^2\right)^{1/2},
\qquad
\varepsilon_W(x,y)=
\frac{\left(\sum_{i\in I_W}|e_i|^2\right)^{1/2}}
{\left(\sum_{i\in I_W}|y_i|^2\right)^{1/2}}.
\]

The maximum absolute error is \(\max_{i\in I_W}|e_i|\). The time window RMS uses the arithmetic mean of squared errors over samples. The rank statistics average across cases, which is a separate operation. The nominal late windows are \([75,150]/\omega_c\) for generic and thermal comparisons and approximately \([100,200]\) fs for BChl. Use the actual stored window endpoints for the robustness data.

| Data family | Full score | Late statistic |
| --- | --- | --- |
| Generic or BChl population | Test minus selected comparison trajectory, normalized by the full reference norm | Relative \(L_2\) on the late samples, normalized by the reference norm on the same samples |
| October finite-temperature comparison | Full-window population relative \(L_2\) | Absolute population RMS for \(t\omega_c\ge75\) |
| November robustness | One row per observable and logged window | Relative \(L_2\), absolute RMS, and maximum absolute difference reported separately. A reference whose RMS is at or below \(10^{-12}\) leaves the relative error missing. |

The robustness coherence column is \(|\rho_{12}|\), and purity is \(\operatorname{Tr}\rho^2\). For paired checks, the second curve named in the row is the denominator for the relative metric. Time-window differences between two finite-bath curves remain paired differences.

## Bath-correlation function

For a zero-temperature spectral model with the documented cutoff interval,

\[
C_{\rm ref}(t)=\int_{\omega_{\rm lo}}^{\omega_{\rm hi}}J(\omega)e^{-i\omega t}\,d\omega,
\qquad C_{\rm disc}(t)=\sum_k\lambda_k^2 e^{-i\omega_k t}.
\]

Use the benchmark frequency and time grids recorded with each batch. For the generic benchmarks, evaluate the continuum sum \(C_{\rm ref}(t_i)=\sum_j J(\omega_j)\Delta_j e^{-i\omega_j t_i}\) on the logarithmic frequency grid, where \(\Delta_j\) are the composite trapezoid weights. The routine `LMD_CORE/common/lmd_composite_weights.m` computes these weights. The discrete weights in mode tables are the squared coupling amplitudes in the convention of [LMD_CORE](../LMD_CORE/README.md). Complex differences are evaluated with the complex magnitude. The BCF metrics use these three distinct normalizations.

1. Full-window relative \(L_2\) is \(\|C_{\rm disc}-C_{\rm ref}\|_{\rm full}/\|C_{\rm ref}\|_{\rm full}\).
2. The late-window gate is \(\operatorname{RMS}_{\rm late}(C_{\rm disc}-C_{\rm ref})/\operatorname{RMS}_{\rm full}(C_{\rm ref})\). A separate budget diagnostic gives the late residual norm over the full reference norm.
3. The maximum absolute error is scaled by \(|C_{\rm ref}(0)|\).

For the generic zero-temperature rank comparison, join the original BCF tables for 20/40/80 and 60 modes and use the canonical ID summary for ID. That calculation uses the archived full-window scores on their recorded grids.
