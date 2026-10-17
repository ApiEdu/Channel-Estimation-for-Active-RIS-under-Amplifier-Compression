"""Figs. 1, 5 and 12, drawn from e_fig1.json, e_fig5.json, e_fig12.json
(numbers produced by matlab/fig1_nmse_floor.m, fig4_gain_mismatch.m, fig7_saleh.m)."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
RED, BLUE, GREEN, ORANGE, PURPLE = "#c00000", "#1f4e79", "#2e7d32", "#e69f00", "#8b6a8b"
M = lambda v: ("0" if round(v) == 0 else f"{v:g}".replace("-", "−"))     # no "-0" labels

# ------------------------------------------------------------------ Fig. 1
d = json.load(open("e_fig1.json")); x = np.array(d["PT_dBm"]); ibo = np.array(d["IBO_dB"])
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.plot(x, d["linear_dB"], "-", c=BLUE, lw=2, label="Linear model (baseline)")
ax.plot(x, d["conv_dB"], "s-", c=RED, lw=1.8, ms=5, label="Conventional LS")
ax.plot(x, d["comp_dB"], "^--", c=GREEN, lw=1.8, ms=5, label="Compensated LS")
ax.plot(x, d["blmmse_dB"], "o-.", c=ORANGE, lw=1.8, ms=5, label="Bussgang-LMMSE")
i = int(np.argmin(d["conv_dB"]))
ax.axvline(x[i], c="0.4", ls=":", lw=1.2, label=f"Conv. LS minimum (IBO = {ibo[i]:.0f} dB)")
ax.set_xlim(-10, 45); ax.set_ylim(-50, 14)
ax.set_xlabel("Pilot transmit power $P_\\mathrm{T}$  [dBm]"); ax.set_ylabel("NMSE  [dB]")
ax.grid(ls=":", alpha=.6); ax.legend(loc="lower left", fontsize=9, framealpha=.95)
top = ax.twiny(); top.set_xlim(ax.get_xlim())
ticks = [20, 10, 0, -10, -20]
pos = np.interp(ticks, ibo[::-1], x[::-1])               # IBO decreases monotonically with P_T
top.set_xticks(pos); top.set_xticklabels([M(t) for t in ticks])
top.set_xlabel("Input back-off  IBO $=A_\\mathrm{sat}^2/(p^2\\sigma_u^2)$  [dB]")
fig.tight_layout(); fig.savefig("../fig_nmse_floor_ibo.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Fig. 5
d = json.load(open("e_fig5.json")); mm = np.array(d["mismatch_dB"])
fig, ax = plt.subplots(figsize=(7.2, 5.0))
ax.axvspan(3, 10.6, color="#f4dcdc", alpha=.7, lw=0)
ax.axvline(0, c="0.5", ls=":", lw=1.2)
ax.plot(mm, d["conv_dB"], "s-", c=RED, lw=1.8, ms=5, label="Conventional LS (no compensation)")
ax.plot(mm, d["fixed_dB"], "o-", c=ORANGE, lw=2, ms=5.5, label="Compensated, calibrated at nominal $A_\\mathrm{sat}$")
ax.plot(mm, d["oracle_dB"], "^--", c=GREEN, lw=1.8, ms=5, label="Compensated, oracle $A_\\mathrm{sat}$")
ax.text(6.8, -14.6, "assumed MORE compression than real:\ncompensation worse than none", ha="center", va="center", color="#8b0000", fontsize=9)
ax.text(-9.8, 5.4, "assumed LESS compression than real:\nestimate shrunk, benign", ha="left", va="center", color="0.25", fontsize=9)
ax.set_xlim(-10.6, 10.6); ax.set_ylim(-17.5, 7.2)
ax.set_xlabel("$A_\\mathrm{sat}$ mismatch: true $-$ assumed  [dB]"); ax.set_ylabel("NMSE at $P_\\mathrm{T}=30$ dBm  [dB]")
ax.grid(ls=":", alpha=.6); ax.legend(loc="lower left", fontsize=9, framealpha=.95)
fig.tight_layout(); fig.savefig("../fig_gain_mismatch.png", dpi=300); plt.close(fig)

# ------------------------------------------------------------------ Fig. 12
d = json.load(open("e_fig12.json")); xi = np.array(d["IBO_dB"]).round(0)
fig, ax = plt.subplots(figsize=(7.0, 5.0))
ax.axvspan(-15, -22.5, color="#f4dcdc", alpha=.7, lw=0)
ax.plot(xi, d["conv_dB"], "s-", c=RED, lw=2, ms=6, label="Conventional LS")
ax.plot(xi, d["abs_kappa_dB"], "d-.", c=PURPLE, lw=2, ms=6, label="Compensated, $|\\kappa|$ only")
ax.plot(xi, d["complex_kappa_dB"], "o-", c=GREEN, lw=2.4, ms=7, label="Compensated, complex $\\kappa$")
k = int(np.argmin(abs(xi)))
y1, y2 = d["abs_kappa_dB"][k], d["complex_kappa_dB"][k]
ax.annotate("", xy=(xi[k], y2 + .25), xytext=(xi[k], y1 - .25), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.5, shrinkA=0, shrinkB=0))
ax.text(xi[k] + .0, y2 - 1.5, f"{abs(y1 - y2):.1f} dB from\nphase correction", ha="center", va="top", color=GREEN, fontsize=9)
ax.text(-18.6, -6.0, "compensation\nno longer helps", ha="center", va="center", color="#8b0000", fontsize=9)
ax.set_xlim(11.5, -22.5); ax.set_ylim(-16, 4.8); ax.set_xticks([10, 5, 0, -5, -10, -15, -20]); ax.set_xticklabels([M(t) for t in [10, 5, 0, -5, -10, -15, -20]])
ax.set_xlabel("Input back-off  [dB]"); ax.set_ylabel("NMSE  [dB]")
ax.grid(ls=":", alpha=.6); ax.legend(loc="upper left", fontsize=9, framealpha=.95)
fig.tight_layout(); fig.savefig("../fig_saleh.png", dpi=300); plt.close(fig)
print("done")

# ------------------------------------------------------------------ Fig. 3
d = json.load(open("e_fig3.json")); K = np.array(d["K"])
fig, ax = plt.subplots(figsize=(7.0, 4.8))
ax.plot(K, d["linear_dB"], "o-", c=BLUE, lw=2, ms=6, label="Linear model")
ax.plot(K, d["comp_dB"], "^--", c=GREEN, lw=2, ms=6, label="Compensated LS (with compression)")
ax.plot(K, d["conv_dB"], "s-", c=RED, lw=2, ms=6, label="Conventional LS (with compression)")
for key, c in (("linear_dB", BLUE), ("comp_dB", GREEN)):
    y0, y1 = d[key][0], d[key][-1]
    ax.plot([K[0], K[-1]], [y0, y0], ls=":", c=c, lw=1.3)              # level at K = 1
    ax.annotate("", xy=(K[-1], y1), xytext=(K[-1], y0),
                arrowprops=dict(arrowstyle="<->", color=c, lw=1.6, shrinkA=0, shrinkB=0))
    ax.text(K[-1] * 1.12, (y0 + y1) / 2, f"{abs(y0 - y1):.1f} dB", color=c, fontsize=11, fontweight="bold", ha="left", va="center")
ax.set_xscale("log", base=2); ax.minorticks_off(); ax.set_xticks(K); ax.set_xticklabels([str(k) for k in K])
ax.set_xlim(0.8, 125); ax.set_ylim(-40, 0)
ax.set_xlabel("Number of BS antennas $K$"); ax.set_ylabel("NMSE at $P_\\mathrm{T}=30$ dBm  [dB]")
ax.grid(ls=":", alpha=.6); ax.legend(loc="lower left", fontsize=9, framealpha=.95)
fig.tight_layout(); fig.savefig("../fig_antenna_scaling.png", dpi=300); plt.close(fig)
print("fig3 done")
