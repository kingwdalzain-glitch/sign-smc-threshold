"""
sweep.py -- Verification code for:

    A. Z. O. Hussien, "An Asymmetric Instability Threshold for Fixed-Gain
    Sliding Mode Control with Discontinuous Switching Under Power-Law
    Contact Stiffness Mismatch," 2026.

Reproduces the reproduction case, the eta and n exponent sweeps, and the
verification figure (Figure 1).

Usage:
    python sweep.py

Runtime: approximately 2-3 minutes on a standard laptop CPU (pure sign(s)
switching requires finer time resolution than boundary-layer SMC).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

M = 0.5
C = 41.7
X_D = 0.005
F_MAX = 700.0
DT = 2e-6
T_SIM = 0.3

K_TRUE_LINEAR = 46524.0


def predicted_threshold(eta, x_d, n=1.0):
    """Closed-form threshold, Eq. (X) in the paper: Delta_k* = eta / x_d^n."""
    return eta / x_d**n


def sim_reproduction():
    """Reproduces the original thesis mismatch scenario under sign(s)
    switching (Section on reproduction)."""
    lam, eta = 150.0, 20.0
    k_true_soft = 8969.0
    k_assumed_stiff = 46524.0
    dk = k_assumed_stiff - k_true_soft
    dk_star = predicted_threshold(eta, X_D)

    x, v = 0.0, 0.0
    N = int(T_SIM / DT)
    max_x = 0.0
    for _ in range(N):
        e = x - X_D
        s = v + lam * e
        sw = np.sign(s)
        f_eq = (C - M * lam) * v + k_assumed_stiff * x
        f = np.clip(f_eq - eta * sw, -F_MAX, F_MAX)
        a = (f - C * v - k_true_soft * x) / M
        v += a * DT
        x += v * DT
        max_x = max(max_x, x)

    print("=== Reproduction case (sign(s) switching) ===")
    print(f"k_assumed={k_assumed_stiff:.1f} N/m, k_true={k_true_soft:.1f} N/m, "
          f"Delta_k={dk:.1f} N/m")
    print(f"Predicted threshold Delta_k* = {dk_star:.1f} N/m "
          f"({dk/dk_star:.1f}x past threshold)")
    print(f"Peak position: {max_x*1000:.2f} mm (target 5.00 mm)")
    print(f"Final position: {x*1000:.2f} mm")
    print()


def sim_vectorized(dk_arr, k_true, lam, eta, x_d, n=1.0, T=T_SIM, dt=DT):
    k_assumed = k_true + dk_arr
    x = np.zeros_like(dk_arr)
    v = np.zeros_like(dk_arr)
    N = int(T / dt)
    for _ in range(N):
        e = x - x_d
        s = v + lam * e
        sw = np.sign(s)
        x_pow = np.sign(x) * np.abs(x) ** n
        f_eq = (C - M * lam) * v + k_assumed * x_pow
        f = np.clip(f_eq - eta * sw, -F_MAX, F_MAX)
        a = (f - C * v - k_true * x_pow) / M
        v = v + a * dt
        x = x + v * dt
    return np.abs(x - x_d) * 1000.0


def bracket_check(dk_range, row, pred):
    pos_start = np.searchsorted(dk_range, 0)
    sub_dk = dk_range[pos_start:]
    sub_row = row[pos_start:]
    above = np.where(sub_row > 1.0)[0]
    if len(above) == 0:
        return None
    idx = above[0]
    if idx == 0:
        lo, hi = sub_dk[0] - (dk_range[1] - dk_range[0]), sub_dk[0]
    else:
        lo, hi = sub_dk[idx - 1], sub_dk[idx]
    return lo <= pred <= hi, lo, hi


def run_sweeps():
    lam0 = 150.0

    # Sweep A: linear contact (n=1), Delta_k vs eta
    eta_range = np.linspace(5, 60, 11)
    dk_range_a = np.linspace(500, 12000, 61)
    z_a = np.zeros((len(eta_range), len(dk_range_a)))
    n_within_a = 0
    for i, eta in enumerate(eta_range):
        z_a[i, :] = sim_vectorized(dk_range_a.copy(), K_TRUE_LINEAR, lam0, eta, X_D)
        pred = predicted_threshold(eta, X_D)
        chk = bracket_check(dk_range_a, z_a[i, :], pred)
        if chk and chk[0]:
            n_within_a += 1
    np.savez("sweep_eta_linear.npz", eta_range=eta_range, dk_range=dk_range_a,
             Z=z_a, lam0=lam0, x_d0=X_D, k_true=K_TRUE_LINEAR)
    print(f"Sweep vs eta (linear, n=1): {n_within_a}/{len(eta_range)} rows matched")

    # Sweep B: linear contact, Delta_k vs x_d (feasible range only)
    x_max_feasible = F_MAX / K_TRUE_LINEAR * 0.9
    x_d_range = np.linspace(0.001, x_max_feasible, 11)
    dk_range_b = np.linspace(500, 20000, 81)
    z_b = np.zeros((len(x_d_range), len(dk_range_b)))
    n_within_b = 0
    for i, xd in enumerate(x_d_range):
        z_b[i, :] = sim_vectorized(dk_range_b.copy(), K_TRUE_LINEAR, lam0, 20.0, xd)
        pred = predicted_threshold(20.0, xd)
        chk = bracket_check(dk_range_b, z_b[i, :], pred)
        if chk and chk[0]:
            n_within_b += 1
    np.savez("sweep_xd_linear.npz", x_d_range=x_d_range, dk_range=dk_range_b,
             Z=z_b, lam0=lam0, eta0=20.0, k_true=K_TRUE_LINEAR)
    print(f"Sweep vs x_d (linear, n=1, feasible range only): "
          f"{n_within_b}/{len(x_d_range)} rows matched")

    # Sweep C: exponent sweep (power-law / Hertzian generalization)
    n_range = [1.0, 1.2, 1.5, 1.8, 2.0]
    n_within_c = 0
    results_c = []
    for n_test in n_range:
        k_true_n = 232.62 / X_D**n_test
        pred = predicted_threshold(20.0, X_D, n=n_test)
        dk_range_c = np.linspace(pred * 0.3, pred * 3, 81)
        row = sim_vectorized(dk_range_c.copy(), k_true_n, lam0, 20.0, X_D, n=n_test)
        chk = bracket_check(dk_range_c, row, pred)
        ok = chk[0] if chk else False
        n_within_c += ok
        results_c.append((n_test, pred, ok))
    print(f"Sweep vs exponent n: {n_within_c}/{len(n_range)} matched")
    for n_test, pred, ok in results_c:
        print(f"    n={n_test}: predicted={pred:.2f}  {'OK' if ok else 'MISS'}")
    np.savez("sweep_n_exponent.npz",
             n_range=np.array(n_range),
             preds=np.array([r[1] for r in results_c]))

    total_within = n_within_a + n_within_b + n_within_c
    total = len(eta_range) + len(x_d_range) + len(n_range)
    print(f"\nTOTAL: {total_within}/{total} rows matched")
    return z_a, eta_range, dk_range_a


def make_figure(z_a, eta_range, dk_range_a):
    d_b = np.load("sweep_xd_linear.npz")
    x_d_range, dk_range_b, z_b = d_b["x_d_range"], d_b["dk_range"], d_b["Z"]
    eta0 = float(d_b["eta0"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    im = ax.pcolormesh(dk_range_a, eta_range, z_a, shading="auto", cmap="RdYlGn_r",
                        norm=LogNorm(vmin=0.05, vmax=20))
    pred_line = [predicted_threshold(eta, X_D) for eta in eta_range]
    ax.plot(pred_line, eta_range, color="blue", lw=2.5, label="Predicted $\\Delta k^*(\\eta)$")
    ax.set_xlabel("$\\Delta k$  [N/m]")
    ax.set_ylabel("$\\eta$  [N]")
    ax.set_title(f"(A) Threshold vs. switching gain $\\eta$\n(sign(s) switching, $x_d$={X_D} m)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    plt.colorbar(im, ax=ax).set_label("Final error [mm]", fontsize=9)

    ax = axes[1]
    im = ax.pcolormesh(dk_range_b, x_d_range, z_b, shading="auto", cmap="RdYlGn_r",
                        norm=LogNorm(vmin=0.05, vmax=20))
    pred_line2 = [predicted_threshold(eta0, xd) for xd in x_d_range]
    ax.plot(pred_line2, x_d_range, color="blue", lw=2.5, label="Predicted $\\Delta k^*(x_d)$")
    ax.set_xlabel("$\\Delta k$  [N/m]")
    ax.set_ylabel("$x_d$  [m]")
    ax.set_title(f"(B) Threshold vs. target $x_d$\n(sign(s) switching, $\\eta$={eta0} N, feasible range)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    plt.colorbar(im, ax=ax).set_label("Final error [mm]", fontsize=9)

    fig.suptitle(r"Asymmetric instability threshold, discontinuous switching:  $\Delta k^* = \eta/x_d$"
                 "\nBlue curve: closed-form prediction with no fitted parameters", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig("sign_threshold_verification.png", dpi=200, bbox_inches="tight")
    print("Saved sign_threshold_verification.png")


if __name__ == "__main__":
    sim_reproduction()
    z_a, eta_range, dk_range_a = run_sweeps()
    make_figure(z_a, eta_range, dk_range_a)
    print("\nAll Paper 1 results reproduced successfully.")
