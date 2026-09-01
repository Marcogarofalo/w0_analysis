import matplotlib.pyplot as plt

# --- Configurazione Stile e Font (Stile pulito, scritte grandi) ---
plt.rcParams.update({
    "text.usetex": False,
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "font.family": "serif",
    "font.serif": ["cmr10"],
    
    "font.size": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    
    # Tick interni e griglia
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "axes.grid": True,
    "axes.grid.which": "major",
    "grid.color": "gray",
    "grid.linestyle": "--",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.5,
})

# --- Dati forniti ---
# groups = ["TOV\n(Group A)", "Cyprus\n(Group C)"   , "Bern\n(Group B)"]
groups = ["(Group A)", "(Group B)", "(Group C)"]
y_values = [609.66, 612.96, 609.54]
y_errors = [5.60, 6.69, 5.12]

final =570.9974+38.64
# final=570.9974 (4.7982) (2.8887) [5.6007]
final_err=5.6
# --- Configurazione Colori e Simboli Personalizzati ---
# Usiamo i tuoi colori (Arancione, Blu, Rosso/Amaranto) e simboli diversi (Cerchio, Quadrato, Triangolo)
colors = ["#f58231", "#4363d8", "#e6194B"] 
markers = ["o", "s", "^"]

# --- Creazione del Grafico ---
fig, ax = plt.subplots(figsize=(8, 6))

# --- Orizontal Band for Final Value ---
# Shaded area between (final - final_err) and (final + final_err)
ax.axhspan(
    final - final_err, 
    final + final_err, 
    color="gray", 
    alpha=0.2, 
    label="Final Error Band"
)

# Plottiamo ogni punto singolarmente per assegnare colore e marker differenti
for i in range(len(groups)):
    ax.errorbar(
        groups[i], 
        y_values[i], 
        yerr=y_errors[i], 
        fmt=markers[i], 
        color=colors[i],
        # ecolor="black",     # Barre d'errore nere per contrasto ottimale
        capsize=5,          
        markersize=10,      # Leggermente più grande per far risaltare le forme diverse
        linewidth=1.5,
        linestyle="none"
    )

# --- Label dell'Asse Y Richiesta (Formattata in LaTeX) ---
ax.set_ylabel(r"$a_{\mu}^{\mathrm{HVP,r}}(I=1,L_{\mathrm{ref}}) \times 10^{10}$")
ax.set_title("Relative Unblinding ($I=1$)", pad=15)

# Contorno del grafico nero e pulito
ax.patch.set_facecolor('white')
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)

# Ottimizzazione spazi per il salvataggio
plt.tight_layout()

# Salvataggio in PDF
plt.savefig("unblinding_plot_colored.pdf", format="pdf")
# plt.show()
