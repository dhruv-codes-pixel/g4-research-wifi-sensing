"""
M3 Member E: statistical utilities for Doppler experiment results.

Functions:
- bootstrap_ci: percentile bootstrap confidence interval (default 10,000 resamples)
- paired_wilcoxon: paired two-sided Wilcoxon signed-rank test
- summary_table: median and bootstrap CI for configuration x method arrays
- wilcoxon_vs_baseline: paired tests against the no-correction baseline

Expected trial-level array shape: (configurations, trials, methods).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

def bootstrap_ci(data, n_resamples=10_000, confidence=0.95,
                 statistic=np.median, seed=2026):
    """Return estimate, lower percentile CI, upper percentile CI."""
    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan
    estimate = float(statistic(x))
    rng = np.random.default_rng(seed)
    boot = np.empty(n_resamples, dtype=float)
    batch_size = 500
    start = 0
    while start < n_resamples:
        count = min(batch_size, n_resamples - start)
        indices = rng.integers(0, x.size, size=(count, x.size))
        boot[start:start + count] = np.median(x[indices], axis=1) if statistic is np.median else [
            statistic(x[idx]) for idx in indices
        ]
        start += count
    alpha = 1.0 - confidence
    lower, upper = np.quantile(boot, [alpha / 2, 1 - alpha / 2])
    return estimate, float(lower), float(upper)

def paired_wilcoxon(method_values, baseline_values):
    """Paired two-sided Wilcoxon signed-rank test; arrays must be trial-aligned."""
    a = np.asarray(method_values, dtype=float)
    b = np.asarray(baseline_values, dtype=float)
    if a.shape != b.shape:
        raise ValueError("Paired arrays must have identical shapes.")
    valid = np.isfinite(a) & np.isfinite(b)
    a, b = a[valid], b[valid]
    if a.size < 2:
        raise ValueError("At least two valid paired observations are required.")
    result = wilcoxon(a, b, alternative="two-sided", zero_method="wilcox", method="auto")
    return float(result.statistic), float(result.pvalue)

def summary_table(values, method_names, n_resamples=10_000, seed=2026):
    """Summarize values shaped (configurations, trials, methods)."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 3 or x.shape[2] != len(method_names):
        raise ValueError("Expected shape (configurations, trials, methods).")
    rows = []
    for ci in range(x.shape[0]):
        for mi, method in enumerate(method_names):
            med, lo, hi = bootstrap_ci(
                x[ci, :, mi], n_resamples=n_resamples, seed=seed + ci * 100 + mi
            )
            rows.append({
                "configuration_index": ci,
                "method": method,
                "n": int(np.isfinite(x[ci, :, mi]).sum()),
                "median": med,
                "ci95_lower": lo,
                "ci95_upper": hi,
            })
    return pd.DataFrame(rows)

def wilcoxon_vs_baseline(values, method_names):
    """Test each correction method against method_names[0] within each configuration."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 3 or x.shape[2] != len(method_names):
        raise ValueError("Expected shape (configurations, trials, methods).")
    rows = []
    for ci in range(x.shape[0]):
        for mi in range(1, len(method_names)):
            stat, p = paired_wilcoxon(x[ci, :, mi], x[ci, :, 0])
            rows.append({
                "configuration_index": ci,
                "comparison": f"{method_names[mi]} vs {method_names[0]}",
                "wilcoxon_statistic": stat,
                "p_value": p,
            })
    return pd.DataFrame(rows)

def plot_summary(table, metric_label, output_path):
    """Plot median and 95% bootstrap CI for each method across configurations."""
    methods = list(table["method"].drop_duplicates())
    configs = list(table["configuration_index"].drop_duplicates())
    x = np.arange(len(configs))
    width = 0.8 / max(1, len(methods))
    fig, ax = plt.subplots(figsize=(11, 5))
    for mi, method in enumerate(methods):
        sub = table[table["method"] == method].sort_values("configuration_index")
        med = sub["median"].to_numpy()
        lo = sub["ci95_lower"].to_numpy()
        hi = sub["ci95_upper"].to_numpy()
        xpos = x + (mi - (len(methods) - 1) / 2) * width
        ax.bar(xpos, med, width, label=method)
        ax.errorbar(xpos, med, yerr=[med - lo, hi - med],
                    fmt="none", capsize=3, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([str(c + 1) for c in configs])
    ax.set_xlabel("Configuration index")
    ax.set_ylabel(metric_label)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
