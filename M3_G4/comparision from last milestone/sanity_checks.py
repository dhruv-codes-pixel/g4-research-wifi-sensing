"""Sanity checks: compare estimated parameters with the TRUE values known to the simulator.
Run: python sanity_checks.py"""
import numpy as np, csi_sim as c

def circ_std_offset_removed(err):
    z = np.exp(1j * err); off = np.angle(z.mean())
    return np.sqrt(np.mean(np.angle(np.exp(1j * (err - off))) ** 2))

rng = np.random.default_rng(7)
T = 50
print("== Check 1: timing / phase estimation error (K=64, gamma=0.9, ideal gain) ==")
res = {m: [[], []] for m in ("Eq21", "Eq22", "Algo4+Eq25")}
for _ in range(T):
    s = c.simulate(K=64, gamma=0.9, rng=rng); hb = s["ht"] / s["g"][:, None]; f = s["f"]
    for name, est in (("Eq21", c.phase_eq21(hb, f)), ("Eq22", c.phase_eq22(hb)), ("Algo4+Eq25", c.phase_algo4_eq25(hb, f))):
        t, ps = est
        terr = (t - s["tau"]); terr = terr - np.median(terr)           # timing is only defined up to a constant offset
        res[name][0].append(np.sqrt(np.mean(terr ** 2)) * 1e9)
        res[name][1].append(circ_std_offset_removed(ps - s["psi"]))
for k, (a, b) in res.items():
    print(f"  {k:12s} timing RMSE = {np.median(a):6.2f} ns   phase RMSE = {np.degrees(np.median(b)):6.2f} deg")

print("\n== Check 2: gain tracking error in dB vs true applied gain (K=64, gamma=0.9) ==")
rng = np.random.default_rng(8); out = {m: [] for m in ("Eq6", "Algo1", "Algo3")}
for _ in range(T):
    s = c.simulate(K=64, gamma=0.9, rng=rng); true_db = 20 * np.log10(s["g"])
    for name, gh in (("Eq6", c.gain_eq6(s["ht"])), ("Algo1", c.gain_algo1(s["ht"])), ("Algo3", c.gain_algo3(s["ht"]))):
        e = 20 * np.log10(gh) - true_db; e -= e.mean()
        out[name].append(np.sqrt(np.mean(e ** 2)))
for k, v in out.items(): print(f"  {k:6s} RMS gain error = {np.median(v):.3f} dB")
print("  (AGC steps are 0.5 dB, true power fluctuation of the channel is about 0.5 dB at K=64)")
