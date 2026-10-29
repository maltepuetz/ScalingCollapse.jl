# Test data

| File | Source | License | Transformation |
|---|---|---|---|
| `IsingMonteCarlo.h5` | 2D Ising model Monte Carlo data (Binder cumulant, magnetization, susceptibility). TODO: document the source. | | |
| `qubitloss/qubitloss_Ic.h5` | M. Pütz et al., "Flow to Nishimori universality in weakly monitored quantum circuits with qubit loss", [arXiv:2505.22720](https://arxiv.org/abs/2505.22720); data from [Zenodo, doi:10.5281/zenodo.15530142](https://doi.org/10.5281/zenodo.15530142) | [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) | Raw `expect_ZZs` samples reduced to the coherent information `I_c` (mean ± standard error) per angle `theta`, system size `L` and parameter `r` with [`test/generators/qubitloss_reduce.jl`](../generators/qubitloss_reduce.jl) |

## `qubitloss/qubitloss_Ic.h5` layout

- `Ls`: system sizes of the phase diagram sweep (8 to 256)
- `thetas`: 15 angles from 0 to π/2
- `theta{i}/L{L}/r`, `Ic`, `Ic_err`, `nsamples`: one curve per angle index `i` and size `L`;
  `theta1` additionally contains `L512`
