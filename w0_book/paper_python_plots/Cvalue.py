import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import markers, rc
import numpy as np


plt.rcParams.update({
    "text.usetex": False,         # Zero crash su Windows
    "mathtext.fontset": "cm",     # Usa il motore matematico Computer Modern
    "axes.formatter.use_mathtext": True, # Tics label in font matematico
    # Associa la famiglia generica 'serif' al font interno corretto di LaTeX
    "font.family": "serif",
    "font.serif": ["cmr10"],      # cmr10 = Computer Modern Roman
    
    # Se preferisci la versione SANS-SERIF di LaTeX, commenta le 2 righe sopra e usa queste:
    # "font.family": "sans-serif",
    # "font.sans-serif": ["cmss10"], # cmss10 = Computer Modern Sans Serif
    # --- Configurazione FONTSIZE (Personalizzabile) ---
    "font.size": 18,          # Dimensione globale di base (se non sovrascritta sotto)
    
    "axes.labelsize": 18,     # Dimensione per i titoli degli assi (X e Y)
    "axes.titlesize": 18,     # Dimensione per il titolo del grafico (ax.set_title)
    
    "xtick.labelsize": 14,    # Dimensione dei numeri sull'asse X
    "ytick.labelsize": 14,    # Dimensione dei numeri sull'asse Y

    # --- Configurazione TICK MINORS ---
    "xtick.minor.visible": True,  # Attiva i tick minori sull'asse X
    "ytick.minor.visible": True,  # Attiva i tick minori sull'asse Y

    # --- Configurazione TICK INTERNI ---
    "xtick.direction": "in",  # Forzza i tick dell'asse X verso l'interno
    "ytick.direction": "in",  # Forza i tick dell'asse Y verso l'interno        
    
    "legend.fontsize": 14,    # Dimensione del testo dentro la legenda
    "figure.titlesize": 18,   # Dimensione del titolo della figura intera (suptitle)

    "errorbar.capsize":5,
    "lines.markeredgewidth":2.0, 
    # --- Configurazione GRIGLIA AUTOMATICA (Major Ticks) ---
    "axes.grid": True,                   # Attiva la griglia di default su tutti i grafici
    "axes.grid.which": "major",          # Applica solo ai ticks principali (major)
    "grid.color": "gray",                # Colore grigio (puoi usare anche codici esadecimali come "#cccccc")
    "grid.linestyle": "--",              # Stile tratteggiato (dashed)
    "grid.linewidth": 0.5,               # Spessore della linea (width 0.5)
    "grid.alpha": 0.7,                   # Opzionale: trasparenza per non appesantire il grafico (da 0 a 1)
})

# Set up the plot style (equivalent to theme_bw)
# plt.style.use("seaborn-v0_8-whitegrid")  # or 'ggplot' depending on preference
# fig, ax = plt.subplots(figsize=(8, 6))
width, height = 800, 600 
fig, ax = plt.subplots(figsize=(width / 100, height / 100))

# -------------------------------------------------------------------------
# Loop 1: ic in (0, -5)
# -------------------------------------------------------------------------
ic_values = [0, -5]
marker_list = ["o", "s"]  # Renamed to avoid overriding matplotlib.markers module

for idx, ic in enumerate(ic_values):
    path_m = f"/home/garofalo/analysis/flow/data_20//fit_all_beta/data_from_wp25_lin_deriv_mc_la_MDs_C{ic}.txt"

    if os.path.exists(path_m):
        df_m = pd.read_csv(path_m, sep=r"\s+")  # sep=r'\s+' handles space/tab delimiters

        # Calculate x and y based on R logic
        x_val = df_m["a[fm]"] * df_m["a[fm]"]
        y_val = df_m["delta_amuc"] / df_m["amul"]
        label_str = f"WP25 $C={ic}$"

        # Error propagation for 
        y_err =  df_m["ddelta_amuc"] / df_m["amul"] + (df_m["delta_amuc"] * df_m["damul"]) / (df_m["amul"] ** 2)

        # Use errorbar instead of scatter to include x-errors (xerr)
        ax.errorbar(
            x_val, y_val, 
            yerr=y_err, 
            fmt=marker_list[idx], 
            label=label_str,
            linestyle='none'  # Prevents drawing lines connecting the points
        )
    else:
        print(f"Warning: File not found {path_m}")

# -------------------------------------------------------------------------
# Section 2: Data from fpi
# -------------------------------------------------------------------------
path_fpi = "/home/garofalo/analysis/flow/data_20/fit_all_beta/data_from_fpi.txt"
if os.path.exists(path_fpi):
    df_fpi = pd.read_csv(path_fpi, sep=r"\s+")

    # Calculate x and y based on R logic
    x_fpi = df_fpi["a[fm]"] * df_fpi["a[fm]"]
    y_fpi = df_fpi["delta_amuc"] / (df_fpi["amul"] )

    # Error propagation for x_fpi = a^2 -> dx = 2 * a * da
    y_fpi_err = df_fpi["ddelta_amuc"] / (df_fpi["amul"] ) + (df_fpi["delta_amuc"] * df_fpi["damul"]) / ((df_fpi["amul"]) ** 2)

    # Use errorbar for the fpi data as well
    ax.errorbar(
        x_fpi, y_fpi, 
        yerr=y_fpi_err, 
        fmt="^", 
        label="FLAG",
        linestyle='none'
    )
else:
    print(f"Warning: File not found {path_fpi}")

# -------------------------------------------------------------------------
# Plot Styling & Limits (Equivalent to your myplotly settings)
# -------------------------------------------------------------------------
ax.set_xlabel(r"$a^2$", fontsize=12)
ax.set_ylabel(r"$\delta\mu_c /\mu_\ell$", fontsize=12)
ax.set_title("")  # Left blank as in your R snippet

# Set x-axis range (equivalent to xrange = c(0, 0.008))
ax.set_xlim(0, 0.008)

# Optional: Uncomment if you want to enforce the y-limits from R comments
# ax.set_ylim(-0.0002, 0.0006)

# Display legend and show plot
ax.legend(frameon=True, facecolor="white", edgecolor="none")
fpi3reg = "Cvalue_plot"
# plt.tight_layout()
plt.savefig(f"{fpi3reg}.pdf", format="pdf")
plt.close()
