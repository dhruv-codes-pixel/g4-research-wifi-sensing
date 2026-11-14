"""Step 8: write LaTeX table bodies (Report/latex/tables/*.tex) from Report/report_data.json (generated from the result CSVs by 07_report_tables.py)."""
import json, re, os
T = json.load(open("../Report/report_data.json")); out = "../Report/latex/tables/"; os.makedirs(out, exist_ok=True)
def esc(s):
    s = str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")
    s = re.sub(r"(?<![\w])-(?=\d)", "$-$", s); return s.replace("−", "$-$")
def tab(name, header, rows, colspec, groupcol=None):
    L = [r"\begin{tabular}{" + colspec + "}", r"\toprule", " & ".join(header) + r" \\", r"\midrule"]
    prev = None
    for r in rows:
        if groupcol is not None and prev is not None and r[groupcol] != prev: L.append(r"\midrule")
        prev = r[groupcol] if groupcol is not None else None
        L.append(" & ".join(esc(c) for c in r) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]; open(out + name + ".tex", "w").write("\n".join(L) + "\n")
tab("inventory", ["Class", "Date", "Cap.", "Frames", "min", r"$\Delta t$ (ms)", r"Gaps $>0.5$\,s", "Windows", "Use"], T["inv"], "llrrrrrrl")
tab("kpi1", ["Class", "$n$", "A0", "A1 (Eq.~21)", "A1w", "A2 (Eq.~22)", "A3 (Alg.~4)"], T["kpi1"], "lrrrrrr")
# phase spread and C_mean compact
cls = ["E_day", "E_night", "S1", "M4"]; meth = ["A0", "A1", "A1w", "A2", "A3"]
def pivot(col):
    rows = []
    for c in cls:
        rows.append([c] + [next(r[col] for r in T["kpi24"] if r[0] == c and r[1] == m) for m in meth])
    return rows
tab("phase", ["Class"] + meth, pivot(2), "lccccc"); tab("cmean", ["Class"] + meth, pivot(3), "lccccc"); tab("clag", ["Class"] + meth, pivot(4), "lccccc")
tab("rule_a", ["Class", "Cap.", "A1", "A1w", "A2", "A3", "A3$-$A1w", "A3$-$A2", "(a) vs A1w"], T["rule_a"], "llrrrrrrl", 0)
tab("rule_b", ["Class", "$V_{\\rm in}$", "$V_{\\rm oos}$", "mean diff.", "max $|$diff$|$"], T["rule_b"], "lrrrr")
tab("occ", ["Method", "Pair", "$V_{\\rm occ}$", "$V_{\\rm emp}$", r"$\Delta V$", "Separation", "AUC"], T["occ"], "llrrrrr", 0)
tab("gain", ["Class", "Jumps/1000", "Range", "Median jump (dB)", "Std all (dB)", "Std in main state (dB)"], T["gain"], "lrrrrr")
tab("eq21fail", ["Class", "Mean \\%", "Range over captures", "Packets"], T["eq21fail"], "lrrr")
tab("time", ["Method", "min (ms)", "max (ms)"], T["time"], "lrr")
tab("synth", ["Channel", r"$\gamma$", "Floor", "A1", "A1w", "A2", "A3"], T["synth"], "lrrrrrr", 0)
print(sorted(os.listdir(out)))
