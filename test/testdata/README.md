# Test data

| File | Source | License | Transformation |
|---|---|---|---|
| `IsingMonteCarlo.h5` | 2D Ising model Monte Carlo data (Binder cumulant, magnetization, susceptibility). TODO: document the source. | | |
| `qubitloss/qubitloss_Ic.h5` | M. Pütz et al., "Flow to Nishimori universality in weakly monitored quantum circuits with qubit loss", [arXiv:2505.22720](https://arxiv.org/abs/2505.22720); data from [Zenodo, doi:10.5281/zenodo.15530142](https://doi.org/10.5281/zenodo.15530142) | [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) | Raw `expect_ZZs` samples reduced to the coherent information `I_c` (mean ± standard error) per angle `theta`, system size `L` and parameter `r` with [`test/generators/qubitloss_reduce.jl`](../generators/qubitloss_reduce.jl) |
| `exact/ising_strip.h5` | Exact correlation length ξ_L/L of the 2D Ising model on a periodic strip (Kaufman transfer-matrix spectrum), L = 8 to 64, 61 points of T ∈ [2.1, 2.45] | generated | [`test/generators/ising_strip.jl`](../generators/ising_strip.jl); layout `Ls`, `L{L}/{x, y}` |
| `exact/curie_weiss.h5` | Exact finite-N Curie–Weiss Binder ratio and m_rms, N = 64 to 1024, 61 points of T ∈ [0.8, 1.2] | generated | [`test/generators/curie_weiss.jl`](../generators/curie_weiss.jl); layout `Ls`, `binder/L{N}/{x, y}`, `m_rms/L{N}/{x, y}` |
| `synthetic/synthetic.h5` | Synthetic finite-size-scaling data with known parameters (x_c = 0.7, ν = 1.2, p3 = 0.3): nine variants (clean, noisy, uneven grid, corrections to scaling, ...), L = 16 to 128 | generated, seed `Xoshiro(20261029)` | [`test/generators/synthetic.jl`](../generators/synthetic.jl); layout `<variant>/Ls`, `<variant>/L{L}/{x, y[, e]}` |
| `percolation/bond_square.h5` | Newman–Ziff bond percolation on the open L × L square lattice: top-bottom spanning probability and largest-cluster fraction, canonical averages (mean ± standard error over samples) at 81 points per L, uniform in z = (p − 1/2)/(1/2)·L^{0.7} ∈ [−2, 2]; L = 16 to 128 with 10⁴ samples, L = 256 with 2·10³ samples | generated, seed `Xoshiro(20261029)` | [`test/generators/percolation.jl`](../generators/percolation.jl); layout `Ls`, `<obs>/L{L}/{x, y, e}` for `obs` ∈ (`spanning`, `strength`) |
| `autoscale/orderParam_L{8,16,32,64,128,256,512}.dat` | 2D site-percolation order parameter (relative size of the largest cluster) from the example data of autoScale.py: O. Melchert, [github.com/omelchert/autoScale](https://github.com/omelchert/autoScale), commit [`a0e98da2098c2b0561edb2c54b9435f31ff884ec`](https://github.com/omelchert/autoScale/tree/a0e98da2098c2b0561edb2c54b9435f31ff884ec) (`ORDER_PARAMETER/`); see also [arXiv:0910.5403](https://arxiv.org/abs/0910.5403) | [MIT](autoscale/LICENSE) (`LICENSE.md` of the repository) | none (copied verbatim); columns p, y, error, `#` lines are comments |

## `qubitloss/qubitloss_Ic.h5` layout

- `Ls`: system sizes of the phase diagram sweep (8 to 256)
- `thetas`: 15 angles from 0 to π/2
- `theta{i}/L{L}/r`, `Ic`, `Ic_err`, `nsamples`: one curve per angle index `i` and size `L`;
  `theta1` additionally contains `L512`
