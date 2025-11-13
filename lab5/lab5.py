#!/usr/bin/env python3
"""
mc_disk_ln.py
Monte-Carlo estimate of I = ∬_{x^2+y^2<1} ln(1/√(x^2+y^2)) dx dy = π/2

Usage:
  python3 mc_disk_ln.py
Produces:
  - results.csv  (columns: method, n, mean_est, mean_se, mean_abs_err, rel_err)
  - error_vs_n.png
"""
import math, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters: modify if needed
N_LIST = [100, 500, 1000, 5000, 10000, 50000]
REPEATS = 50
SEED = 123456

# exact value
I_exact = math.pi / 2.0

def sample_uniform_disk(n, rng):
    """
    Return array of r values (radius) for n points uniformly in unit disk.
    We actually only need r to compute f = -ln r.
    """
    u = rng.random(n)            # U(0,1)
    r = np.sqrt(u)               # r distribution for uniform disk
    theta = rng.random(n) * 2 * math.pi  # not needed if only r used
    return r

def estimate_once_uniform(n, rng):
    r = sample_uniform_disk(n, rng)
    f = -np.log(r)               # f = -ln r
    mean_f = f.mean()
    var_f = f.var(ddof=1)
    se = math.sqrt(var_f / n)
    I_hat = math.pi * mean_f
    se_I = math.pi * se
    return I_hat, se_I

def run_experiment(n_list, repeats, seed):
    rng_master = np.random.default_rng(seed)
    rows = []
    for n in n_list:
        ests = np.zeros(repeats)
        ses = np.zeros(repeats)
        for r in range(repeats):
            sub_seed = rng_master.integers(0, 2**31-1)
            rng = np.random.default_rng(sub_seed)
            Ihat, seI = estimate_once_uniform(n, rng)
            ests[r] = Ihat
            ses[r] = seI
        mean_est = ests.mean()
        mean_se = ses.mean()
        mean_abs_err = np.mean(np.abs(ests - I_exact))
        rel_err = mean_abs_err / abs(I_exact)
        rows.append({
            "method": "uniform_disk",
            "n": n,
            "mean_est": mean_est,
            "mean_se": mean_se,
            "mean_abs_err": mean_abs_err,
            "rel_err": rel_err
        })
        print(f"n={n}: mean_est={mean_est:.6f}, mean_abs_err={mean_abs_err:.6f}, mean_se={mean_se:.6f}")
    df = pd.DataFrame(rows)
    return df

def plot_results(df, png_out="error_vs_n.png"):
    pivot = df.pivot(index="n", columns="method", values="mean_abs_err")
    plt.figure(figsize=(7,5))
    for col in pivot.columns:
        plt.plot(pivot.index, pivot[col], marker='o', label=col)
    ns = np.array(sorted(pivot.index))
    # reference 1/sqrt(n) scaled
    ref = pivot.max(axis=1).iloc[0] * (math.sqrt(ns[0]) / np.sqrt(ns))
    plt.plot(ns, ref, linestyle='--', color='gray', label='~ 1/sqrt(n)')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel('n (log scale)'); plt.ylabel('mean abs error (log scale)')
    plt.title('MC error vs n (uniform disk sampling)')
    plt.legend(); plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout(); plt.savefig(png_out)
    print(f"Saved plot {png_out}")

def main():
    t0 = time.time()
    print("Exact I =", I_exact)
    df = run_experiment(N_LIST, REPEATS, SEED)
    df.to_csv("results.csv", index=False)
    print("Saved results.csv")
    plot_results(df, "error_vs_n.png")
    print("Done in {:.1f}s".format(time.time()-t0))

if __name__ == "__main__":
    main()
