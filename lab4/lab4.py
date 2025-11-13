#!/usr/bin/env python3
"""
Monte-Carlo estimation of I = ∫_{-∞}^{∞} e^{-x^6} dx

Produces:
 - CSV file with results for different n and two methods (Uniform truncation and Importance Sampling)
 - PNG plot error_vs_n (log-log)
"""
import os
import math
import argparse
from time import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import special
from tqdm import trange

# -----------------------
# Parameters (you can change these)
# -----------------------
N_LIST = [100, 500, 1000, 5000, 10000, 50000]   # values of n to test
REPEATS = 50             # number of independent repetitions for each n
L_TRUNC = 3.0            # truncation limit for uniform method [-L,L]
SIGMA_IS = 1.0           # sigma for Normal proposal in importance sampling
SEED = 123456            # base RNG seed for reproducibility

OUT_CSV = "results.csv"
OUT_PNG = "error_vs_n.png"

# -----------------------
# Exact value
# -----------------------
I_exact = (1.0 / 3.0) * special.gamma(1.0 / 6.0)  # exact: (1/3)*Gamma(1/6)


# -----------------------
# Helpers
# -----------------------
def estimate_uniform_truncation(n, L, rng):
    """
    One Monte-Carlo estimate using uniform sampling on [-L, L].
    Returns (estimate, estimated_se)
    """
    u = rng.uniform(-L, L, size=n)
    fx = np.exp(-np.power(u, 6))
    vals = (2.0 * L) * fx               # importance weight: volume * fx (since density = 1/(2L))
    mean = vals.mean()
    # sample variance (unbiased)
    var = vals.var(ddof=1)
    se = math.sqrt(var / n)
    return mean, se


def estimate_importance_normal(n, sigma, rng):
    """
    One Monte-Carlo importance-sampling estimate using Normal(0, sigma^2) as proposal.
    Returns (estimate, estimated_se)
    """
    x = rng.normal(loc=0.0, scale=sigma, size=n)
    fx = np.exp(-np.power(x, 6))
    # proposal density phi(x)
    phi = (1.0 / (math.sqrt(2.0 * math.pi) * sigma)) * np.exp(- (x * x) / (2.0 * sigma * sigma))
    w = fx / phi
    mean = w.mean()
    var = w.var(ddof=1)
    se = math.sqrt(var / n)
    return mean, se


# -----------------------
# Experiment runner
# -----------------------
def run_experiment(n_list, repeats, L, sigma, seed):
    rng_master = np.random.default_rng(seed)
    rows = []
    total_runs = len(n_list) * repeats * 2
    run_idx = 0
    t0 = time()

    for n in n_list:
        # arrays to collect for this n
        uni_est = np.zeros(repeats)
        uni_se = np.zeros(repeats)
        is_est = np.zeros(repeats)
        is_se = np.zeros(repeats)

        for r in range(repeats):
            # create a sub-RNG for reproducibility per repetition
            subseed = rng_master.integers(0, 2**31 - 1)
            rng = np.random.default_rng(subseed)

            est_u, se_u = estimate_uniform_truncation(n, L, rng)
            uni_est[r] = est_u
            uni_se[r] = se_u
            run_idx += 1

            # importance sampling
            # use different RNG stream
            subseed2 = rng_master.integers(0, 2**31 - 1)
            rng2 = np.random.default_rng(subseed2)
            est_is, se_is = estimate_importance_normal(n, sigma, rng2)
            is_est[r] = est_is
            is_se[r] = se_is
            run_idx += 1

        # compute summary statistics across repeats
        row_u = {
            "method": "uniform_trunc",
            "n": n,
            "mean_est": uni_est.mean(),
            "std_est": uni_est.std(ddof=1),
            "mean_se": uni_se.mean(),
            "mean_abs_err": np.mean(np.abs(uni_est - I_exact)),
            "rel_err": np.mean(np.abs(uni_est - I_exact)) / abs(I_exact)
        }
        rows.append(row_u)

        row_is = {
            "method": "importance_normal",
            "n": n,
            "mean_est": is_est.mean(),
            "std_est": is_est.std(ddof=1),
            "mean_se": is_se.mean(),
            "mean_abs_err": np.mean(np.abs(is_est - I_exact)),
            "rel_err": np.mean(np.abs(is_est - I_exact)) / abs(I_exact)
        }
        rows.append(row_is)

        # progress print
        elapsed = time() - t0
        print(f"n={n}: uniform mean={row_u['mean_est']:.6f}, IS mean={row_is['mean_est']:.6f}, elapsed={elapsed:.1f}s")

    df = pd.DataFrame(rows)
    return df


# -----------------------
# Plotting
# -----------------------
def plot_results(df, out_png):
    # pivot to plot mean_abs_err vs n for each method
    pivot = df.pivot(index="n", columns="method", values="mean_abs_err")
    plt.figure(figsize=(7,5))
    for col in pivot.columns:
        plt.plot(pivot.index, pivot[col], marker='o', label=col)
    # reference 1/sqrt(n) line (scaled)
    ns = np.array(sorted(pivot.index))
    ref = pivot.max(axis=1).iloc[0] * (np.sqrt(ns[0]) / np.sqrt(ns))  # scale to start near first point
    plt.plot(ns, ref, linestyle='--', color='gray', label=r'~ 1/sqrt(n)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('n (log scale)')
    plt.ylabel('mean absolute error (log scale)')
    plt.title('MC error vs n')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved plot to {out_png}")


# -----------------------
# Main
# -----------------------
def main():
    parser = argparse.ArgumentParser(description="Monte Carlo estimation of integral ∫ e^{-x^6} dx")
    parser.add_argument('--out', default=OUT_CSV, help="CSV output file")
    parser.add_argument('--png', default=OUT_PNG, help="PNG plot file")
    parser.add_argument('--repeats', type=int, default=REPEATS, help="repeats per n")
    parser.add_argument('--seed', type=int, default=SEED, help="random seed")
    args = parser.parse_args()

    print("Exact I =", I_exact)
    df = run_experiment(N_LIST, args.repeats, L_TRUNC, SIGMA_IS, args.seed)
    df.to_csv(args.out, index=False)
    print(f"Saved results CSV to {args.out}")
    plot_results(df, args.png)
    print("Done.")

if __name__ == "__main__":
    main()
