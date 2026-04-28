"""Build FIGURE_S3_panel_signal_synthesis.

Single supplementary figure (option (b) for the cross-application synthesis):
plots panel-level Cramer's V vs hierarchy depth for the two applications.

Numerical inputs:
- CT full-VOG outcome association: V at depths 1-4 from the CT full-VOG
  hierarchy results referenced by the combined manuscript (depth-4 V = 0.371
  is reported in the manuscript main text; depths 1-3 are from the
  CT_clearance analysis output table cited by the manuscript).
- AGP self-reported IBD exploratory condition-prioritization: V at depth 1
  and depth 4 (0.24, 0.31) are reported in the main text. Depths 2 and 3
  are not reported in the main text and are shown as a dashed-line
  endpoint-to-endpoint guide rather than imputed.

Note: AGP IBD prioritization is an exploratory phenotype screen, NOT the
Halfvarson direct-transfer / rebuilt-cohort test that supports the gut
portability claim. Halfvarson panel-level signal is reported categorically
(q < 0.05) in the main text rather than as a depth-resolved V curve.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT_DIR = "/sessions/awesome-confident-turing/mnt/dcst-methods/assets/figures"
os.makedirs(OUT_DIR, exist_ok=True)

# CT full-VOG panel-level Cramer's V across depths 1-4
ct_depths = [1, 2, 3, 4]
ct_v      = [0.214, 0.253, 0.357, 0.371]

# AGP self-reported IBD exploratory phenotype prioritization V endpoints
agp_depths_anchors = [1, 4]
agp_v_anchors      = [0.24, 0.31]

fig, ax = plt.subplots(figsize=(6.0, 4.0))

# CT full-VOG
ax.plot(ct_depths, ct_v, marker="o", color="#1f6f8b", linewidth=2,
        label="CT full-VOG outcome panel ($n_{CT}=604$)")
for x, y in zip(ct_depths, ct_v):
    ax.annotate(f"{y:.3f}", (x, y), textcoords="offset points",
                xytext=(6, 6), fontsize=8, color="#1f6f8b")

# AGP IBD prioritization (endpoints only, dashed guide, no imputation)
ax.plot(agp_depths_anchors, agp_v_anchors, marker="s", color="#c0492b",
        linewidth=1.4, linestyle="--",
        label="AGP self-reported IBD condition screen (depth-1 and depth-4 endpoints reported)")
for x, y in zip(agp_depths_anchors, agp_v_anchors):
    ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                xytext=(6, -12), fontsize=8, color="#c0492b")

ax.set_xticks([1, 2, 3, 4])
ax.set_xlabel("Hierarchy depth $\\ell$")
ax.set_ylabel("Panel-level Cramer's $V$")
ax.set_ylim(0.10, 0.42)
ax.set_xlim(0.7, 4.3)
ax.grid(True, axis="y", alpha=0.3, linestyle=":")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", fontsize=8, frameon=False)

ax.text(0.99, 0.02,
        ("Scales and tests differ by design: AGP IBD is an exploratory phenotype\n"
         "screen in the AGP discovery cohort; CT full-VOG is the primary outcome\n"
         "association test in the labeled CT subset. Halfvarson direct-transfer\n"
         "and rebuilt-cohort panel-level signal is reported categorically in the\n"
         "main text (Figure 3) rather than as a depth-resolved $V$ curve."),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
        color="#444444",
        bbox=dict(facecolor="white", edgecolor="#cccccc", boxstyle="round,pad=0.3"))

plt.tight_layout()

out_pdf = os.path.join(OUT_DIR, "FIGURE_S3_panel_signal_synthesis.pdf")
out_png = os.path.join(OUT_DIR, "FIGURE_S3_panel_signal_synthesis.png")
plt.savefig(out_pdf, bbox_inches="tight")
plt.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")
