"""Fig. 13: plotted from e_spread.json (numbers produced by matlab/fig8_element_spread.m)."""
import json, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = json.load(open("e_spread.json")); x = d["sigma_A_dB"]
fig, ax = plt.subplots(figsize=(7, 4.6)); 
ax.plot(x, d["conventional_LS_dB"], "s-", c="#c00000", lw=1.8, label="Conventional LS")
ax.plot(x, d["single_kappa_nominal_dB"], "d-", c="#8b6a8b", lw=1.8, label="Single $\\kappa$, nominal value")
ax.plot(x, d["single_kappa_mean_dB"], "o-", c="#e69f00", lw=2.0, label="Single $\\kappa$, set to the mean")
ax.plot(x, d["per_element_oracle_dB"], "^--", c="#2e7d32", lw=1.8, label="Per-element $\\kappa_n$ (oracle)")
gap = d["single_kappa_nominal_dB"][-1] - d["per_element_oracle_dB"][-1]
ax.annotate("", xy=(4, d["per_element_oracle_dB"][-1]), xytext=(4, d["single_kappa_nominal_dB"][-1]),
            arrowprops=dict(arrowstyle="<->", color="#555", lw=1.6))
ax.text(3.88, (d["per_element_oracle_dB"][-1] + d["single_kappa_nominal_dB"][-1]) / 2,
        f"{abs(gap):.1f} dB", ha="right", va="center", fontsize=11, fontweight="bold", color="#333")
ax.set_xlabel("Per-element $A_{\\mathrm{sat}}$ spread, $\\sigma_A$  [dB]")
ax.set_ylabel("NMSE at $P_T=30$ dBm  [dB]"); ax.grid(ls=":", alpha=.6)
ax.legend(loc="center left", fontsize=9, framealpha=.95)
fig.tight_layout(); fig.savefig("../fig_element_spread.png", dpi=300)
print("gap nominal-oracle at 4 dB: %.2f dB" % abs(gap))
