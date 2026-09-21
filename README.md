# An Asymmetric Instability Threshold for Fixed-Gain Sliding Mode Control with Discontinuous Switching Under Power-Law Contact Stiffness Mismatch

Verification code and results for:

> A. Z. O. Hussien, "An Asymmetric Instability Threshold for Fixed-Gain Sliding Mode Control with Discontinuous Switching Under Power-Law Contact Stiffness Mismatch," 2026. [arXiv link to be added on posting]

## What this is

This repository contains the complete, self-contained simulation code used to derive and verify this paper's central result: a closed-form sliding-existence threshold

```
Delta_k* = eta / x_d^n
```

for a classical discontinuous-switching (sign(s)) sliding mode controller regulating a single-degree-of-freedom plant with power-law contact force `F = k*|x|^n*sign(x)`. This covers both linear contact (n=1) and Hertzian contact (n=3/2, the standard model for compliant and biological tissue contact).

Overestimating the stiffness coefficient past this threshold destroys the sliding condition at the target and drives the system to a saturation-locked false equilibrium. Underestimating it leaves ideal sliding-mode tracking exact.

This threshold is the exact `Phi -> 0` limit of a companion result for boundary-layer (saturation-smoothed) SMC — see [asymmetric-smc-threshold](https://github.com/azoh-controls/asymmetric-smc-threshold) — obtained independently via a different derivation technique (classical equivalent-control / sliding-existence analysis here, versus a boundary-layer fixed-point argument there). The exact agreement between the two is a nontrivial cross-check of both results.

## Reproducing the paper's results

```bash
pip install numpy matplotlib
python sweep.py
```

Runtime: approximately 2-3 minutes (discontinuous switching requires finer time resolution than boundary-layer SMC).

This reproduces, in order:

1. **The original failure case** — a controller tuned for a stiff object and deployed on a much softer one, reproducing the reported 1,941%-scale overshoot behavior.
2. **Three independent verification sweeps**: switching gain `eta`, target displacement `x_d` (restricted to the actuator-feasible range), and contact exponent `n` — 27 total configurations, all matching the closed-form prediction exactly.
3. **Figure 1**, saved as `sign_threshold_verification.png`, exactly as it appears in the paper.

## Files

| File | Description |
|---|---|
| `sweep.py` | All simulation, verification, and plotting code |
| `sweep_eta_linear.npz` | Data for Figure 1, panel A |
| `sweep_xd_linear.npz` | Data for Figure 1, panel B |
| `sweep_n_exponent.npz` | Data for the contact-exponent verification sweep |
| `sign_threshold_verification.png` | Figure 1, as it appears in the paper |

## Verification summary

27/27 tested configurations (11 across switching gain, 11 across target displacement, 5 across contact exponent) match the closed-form prediction exactly, with no exceptions.

## Model

The plant is a single-degree-of-freedom mass-damper with power-law contact,

```
m*x_ddot + c*x_dot + k_true*|x|^n*sign(x) = F
```

controlled by a classical discontinuous sliding mode law with sliding surface `s = x_dot + lambda*(x - x_d)` and control

```
F = sat( (c - m*lambda)*x_dot + k_assumed*|x|^n*sign(x) - eta*sign(s), F_max )
```

See the paper for the full derivation via the classical equivalent-control (Utkin) method.

## Citation

```bibtex
@article{hussien2026sign,
  title   = {An Asymmetric Instability Threshold for Fixed-Gain Sliding Mode Control with Discontinuous Switching Under Power-Law Contact Stiffness Mismatch},
  author  = {Hussien, Ahmed Z. O.},
  year    = {2026},
  note    = {arXiv preprint}
}
```

## License

MIT License — see [LICENSE](LICENSE).
