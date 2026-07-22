# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# # 1. Read data (handles arbitrary whitespace fields automatically)
# df = pd.read_csv("fpi.txt", sep=r"\s+")

# # 2. Setup structural plotting panels
# fig, axes = plt.subplots(1, 2, figsize=(11, 5))
# ensembles = df['en'].unique()

# for ax, en in zip(axes, ensembles):
#     sub = df[df['en'] == en].copy().reset_index()
    
#     # Define horizontal groupings for clean visualization
#     # Left group: finite values side-by-side | Right group: infinite values side-by-side
#     x_finite = [0.8, 1.2]
#     x_infinite = [1.8, 2.2]
    
#     # Plot Finite volume points
#     ax.errorbar(x_finite, sub['value'], yerr=sub['error'], 
#                 fmt='o', color='crimson', capsize=5, elinewidth=1.5, label='Finite L')
    
#     # Plot Infinite volume points
#     ax.errorbar(x_infinite, sub['Linf_value'], yerr=sub['Linf_error'], 
#                 fmt='s', color='navy', capsize=5, elinewidth=1.5, label=r'Extrapolated $L_{\infty}$')
    
#     # Visual markers: Annotate the lattice size next to the finite points
#     for i, txt in enumerate(sub['L']):
#         ax.annotate(f"L={txt}", (x_finite[i], sub['value'][i]), textcoords="offset points", 
#                     xytext=(-25, -3), ha='center', fontsize=9, color='crimson')

#     # Formatting tweaks
#     ax.set_title(f"Ensemble {en}", fontsize=13, fontweight='bold')
#     ax.set_xticks([1.0, 2.0])
#     ax.set_xticklabels(['Finite Sizes', r'Infinite Vol. ($L_{\infty}$)'], fontsize=11)
#     ax.grid(True, linestyle=':', alpha=0.6)
    
# axes[0].set_ylabel(r"$f_{\pi}$ decay constant", fontsize=12)
# axes[0].legend(loc='best')
# plt.tight_layout()
# plt.savefig("fpi_extrapolation_check.png", dpi=300)
# plt.show()
 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
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

# Read your dataset
df = pd.read_csv("fpi_20.txt", sep=r"\s+")

ensembles = df['en'].unique()
fig, axes = plt.subplots(1, len(ensembles), figsize=(12, 6), sharey=False)

# Define explicit color palettes for distinct volumes per ensemble
color_maps = {
    # 'B': {64: '#e66101', 96: '#5e3c99'},
    # 'C': {80: '#d7191c', 112: '#2c7bb6'},
    # 'D': {96: '#fdae61', 128: '#abd9e9'} # Example color palette extension
    'B': {64: '#d7191c', 96: '#2c7bb6'},
    'C': {80: '#d7191c', 112: '#2c7bb6'},
    'D': {96: '#d7191c', 128: '#2c7bb6'} # Example color palette extension
}
marker_maps = {
    'B': {64: 'o', 96: 's'},
    'C': {80: 'o', 112: 's'},
    'D': {96: 'o', 128: 's'} # Example color palette extension
}
markerL_maps = {
    'B': {64: '^', 96:  'x'},
    'C': {80: '^', 112: 'x'},
    'D': {96: '^', 128: 'x'} # Example color palette extension
}

for ax, en in zip(axes, ensembles):
    sub = df[df['en'] == en]
    
    i=0
    for _, row in sub.iterrows():
        L_size = int(row['L'])
        point_color = color_maps[en][L_size]
        marker_style = marker_maps[en][L_size]
        markerL_style = markerL_maps[en][L_size]
        x_finite = [0.9, 1.1]
        x_infinite = [1.9, 2.1]
        # Plot finite-volume point (Circle)
        ax.errorbar(x_finite[i], row['value'], yerr=row['error'], 
                    fmt=marker_style,
                    color=point_color,
                    ms=8, capsize=5, 
                    label=f"L = {L_size}")
        
        # Plot its corresponding infinite-volume point (Square)
        ax.errorbar(x_infinite[i], row['Linf_value'], yerr=row['Linf_error'], 
                    fmt=markerL_style,
                    color=point_color,
                    ms=8, capsize=5
                    )
        
        # Draw explicit linking trajectory line
        ax.plot([x_finite[i], x_infinite[i]], [row['value'], row['Linf_value']], 
                linestyle='--', color=point_color, alpha=0.7, linewidth=1.5)
        i+=1
    # Style configuration
    ax.set_title(f"Ensemble {en}", fontsize=13, fontweight='bold')
    ax.set_xticks([1.0, 2.0])
    ax.set_xticklabels(['Finite $L$', r'Extrapolated to $L_{\infty}$'], fontsize=11)
    ax.set_xlim(0.6, 2.4)
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))  # Format Y-axis to 4 decimal places
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(title="Lattice Volume", loc='lower right')

axes[0].set_ylabel(r"$af_{\pi}^{\rm WTI}$", fontsize=12)
plt.tight_layout()
plt.tight_layout()

# 4. Save directly as a single-page vector graphic PDF
plt.savefig("fpi_FVE_plot.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.close()

print("Success! Your plot has been saved to 'fpi_FVE_plot.pdf'.")
