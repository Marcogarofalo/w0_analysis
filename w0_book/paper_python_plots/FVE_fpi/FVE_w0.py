import matplotlib.pyplot as plt
import pandas as pd

# 1. Read the data file (handling whitespace separation safely)
df = pd.read_csv("w0.txt", sep=r"\s+")

# 2. Get unique ensembles to determine the layout
ensembles = df["en"].unique()
num_plots = len(ensembles)

# 3. Create a side-by-side subplot layout
fig, axes = plt.subplots(
    1, num_plots, figsize=(6 * num_plots, 5), sharex=False, squeeze=False
)
a_fm = [0.079696, 0.056879]
# 4. Loop through each ensemble and plot its data
for idx, en_name in enumerate(ensembles):
    ax = axes[0, idx]
    # Filter data for the specific ensemble
    sub_df = df[df["en"] == en_name].sort_values("L")

    # Calculate scale values
    x_data = sub_df["L"] * a_fm[idx]
    y_data = sub_df["value"]
    y_err = sub_df["error"]

    # Plot values with error bars
    ax.errorbar(
        x_data,
        y_data,
        yerr=y_err,
        fmt="o",
        capsize=5,
        markersize=6,
        label=f"Ensemble {en_name}",
    )

    # Styling and labeling
    ax.set_title(f"Ensemble: {en_name}", fontsize=14)
    ax.set_xlabel("$L$ [fm]", fontsize=12)
    ax.set_ylabel("$w_0/a$ after $m_0$ correction", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    #ax.legend()
    # ax.legend(title="Lattice Volume", loc='lower right')

    # --- ADJUST RANGES HERE ---
    # Enlarge X range by adding a 15% margin on both sides
    x_min, x_max = x_data.min(), x_data.max()
    x_margin = (x_max - x_min) * 0.15 if x_max != x_min else 0.1
    ax.set_xlim(x_min - x_margin, x_max + x_margin)

    # Enlarge Y range considering error bars (15% margin)
    y_min = (y_data - y_err).min()
    y_max = (y_data + y_err).max()
    y_margin = (y_max - y_min) * 0.15 if y_max != y_min else 0.1
    ax.set_ylim(y_min - y_margin, y_max + y_margin)
# --------------------------
# 5. Optimize spacing and render
# plt.tight_layout()
# plt.show()
plt.savefig("w0_FVE_plot.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.close()

