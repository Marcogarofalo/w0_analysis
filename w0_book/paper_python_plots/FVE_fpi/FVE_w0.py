import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.container as mcontainer
from matplotlib.legend_handler import HandlerErrorbar
 
handler_map = {
    mcontainer.ErrorbarContainer: HandlerErrorbar(yerr_size=.7)
}
def legend(ax, *args, **kwargs):
    kwargs.setdefault("handler_map", handler_map)
    return ax.legend(*args, **kwargs)

# Configurazione corretta per usare i font interni di LaTeX ovunque
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
    
    "legend.fontsize": 14,    # Dimensione del testo dentro la legenda
    "figure.titlesize": 18,   # Dimensione del titolo della figura intera (suptitle)

    "legend.markerscale": 1.,       # Ingrandisce i simboli solo dentro la legenda (moltiplicatore)
    "legend.labelspacing": 1.1,      # Aumenta lo spazio verticale tra le righe (default 0.5)
    "legend.borderpad": .8,
    "legend.handletextpad": 1.,


    "errorbar.capsize":5,
    "lines.markeredgewidth":2.0, 
    "lines.markersize": 9.0, 
    # --- Configurazione TICK MINORS ---
    "xtick.minor.visible": True,  # Attiva i tick minori sull'asse X
    "ytick.minor.visible": True,  # Attiva i tick minori sull'asse Y
    # --- Configurazione TICK INTERNI ---
    "xtick.direction": "in",  # Forzza i tick dell'asse X verso l'interno
    "ytick.direction": "in",  # Forza i tick dell'asse Y verso l'interno        
 
    # --- Configurazione GRIGLIA AUTOMATICA (Major Ticks) ---
    "axes.grid": True,                   # Attiva la griglia di default su tutti i grafici
    "axes.grid.which": "major",          # Applica solo ai ticks principali (major)
    "grid.color": "gray",                # Colore grigio (puoi usare anche codici esadecimali come "#cccccc")
    "grid.linestyle": "--",              # Stile tratteggiato (dashed)
    "grid.linewidth": 0.5,               # Spessore della linea (width 0.5)
    "grid.alpha": 0.7,                   # Opzionale: trasparenza per non appesantire il grafico (da 0 a 1)
})


blue = "#4363d8"
orange = "#f58231"
yellow = "#ffe119"
maroon = "#800000"
navy = "#000075"
lavender = "#dcbeff"
red = "#e6194B"


# 1. Read the data file (handling whitespace separation safely)
df = pd.read_csv("w0_20.txt", sep=r"\s+")

# 2. Get unique ensembles to determine the layout
ensembles = df["en"].unique()
num_plots = len(ensembles)

# 3. Create a side-by-side subplot layout
fig, axes = plt.subplots(
    1, num_plots, figsize=(12, 6), sharex=False
)
a_fm = [0.079696, 0.068004, 0.056879]
# 4. Loop through each ensemble and plot its data
for idx, en_name in enumerate(ensembles):
    ax = axes[idx]
    
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
    ax.set_title(f"Ensemble: {en_name}")
    ax.set_xlabel("$L$ [fm]")
    ax.set_ylabel("$w_0/a$ after $m_0$ correction")
    if idx > 0:
        ax.yaxis.set_label_text("")
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
# plt.show()


plt.tight_layout()
plt.savefig("w0_FVE_plot.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.close()

