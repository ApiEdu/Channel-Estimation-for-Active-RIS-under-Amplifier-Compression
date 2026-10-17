"""Fig. 2 from e_training_amp.json (numbers from matlab/fig2_training_amplitude.m)."""
import json, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = json.load(open("e_training_amp.json")); a = d["a"]; lin = d["linear_dB"]; nl = d["compressed_dB"]
im = min(range(len(nl)), key=lambda i: nl[i]); ao, yo = a[im], nl[im]
pen = nl[-1] - yo
fig, ax = plt.subplots(figsize=(7, 4.8))
ax.plot(a, lin, "o-", c="#1f4e79", lw=2, label="Linear model (prior-art design assumption)")
ax.plot(a, nl, "s-", c="#c00000", lw=2, label="With amplifier compression")
ax.plot([ao], [yo], "*", c="#2e7d32", ms=16, zorder=5)
ax.annotate(f"true optimum\n$a={ao:g}$", xy=(ao, yo), xytext=(4.6, -1.2), color="#2e7d32", fontsize=10,
            ha="left", arrowprops=dict(arrowstyle="->", color="#2e7d32", lw=1.4))
# reference level of the optimum and the penalty at a = 10
ax.plot([ao, a[-1]], [yo, yo], ls=":", c="#c00000", lw=1.3)
ax.annotate("", xy=(a[-1], nl[-1]), xytext=(a[-1], yo),
            arrowprops=dict(arrowstyle="<->", color="#c00000", lw=1.6, shrinkA=0, shrinkB=0))
ax.text(a[-1] - 0.15, (yo + nl[-1]) / 2, f"{pen:.1f} dB", color="#c00000", fontsize=11,
        fontweight="bold", ha="right", va="center")
ax.annotate("linear model says\n“use maximum gain”", xy=(a[-1] - 0.05, lin[-1] + 0.4), xytext=(7.0, -10.3),
            color="#1f4e79", fontsize=10, ha="left",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.3))
ax.set_xlim(0.6, 10.8); ax.set_ylim(-21.5, 3.2); ax.set_xticks([1, 2, 4, 6, 8, 10])
ax.set_xlabel("Training amplification $a$")
ax.set_ylabel("Sum-NMSE over $[\\mathbf{h}_\\mathrm{d};\\,\\mathbf{b}]$ at $P_T=30$ dBm  [dB]")
ax.grid(ls=":", alpha=.6); ax.legend(loc="lower left", fontsize=9, framealpha=.95)
fig.tight_layout(); fig.savefig("../fig_training_amplitude.png", dpi=300)
print("optimum a=%g, penalty %.2f dB" % (ao, pen))
