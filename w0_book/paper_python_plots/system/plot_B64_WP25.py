import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib import rc

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

def parse_jackknife_with_errors(filepath):
    """
    Parses central values and jackknife uncertainties from the system data file.
    """
    data = {}
    current_key = None
    buffer_lines = []

    # Regex matching decimals and numbers in scientific notation
    num_pattern = re.compile(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\b[-+]?\d+\b")

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Section header switches
            if "Matrix" in line:
                current_key = "matrix"
                buffer_lines = []
                continue
            elif "Riso:" in line:
                if current_key and buffer_lines: data[current_key] = buffer_lines
                current_key = "Riso"
                buffer_lines = []
                continue
            elif "Rsim:" in line:
                if current_key and buffer_lines: data[current_key] = buffer_lines
                current_key = "Rsim"
                buffer_lines = []
                continue
            elif "solution:" in line:
                if current_key and buffer_lines: data[current_key] = buffer_lines
                current_key = "solution"
                buffer_lines = []
                continue
            
            matches = num_pattern.findall(line)
            if matches:
                buffer_lines.append([float(x) for x in matches])
        
        if current_key and buffer_lines:
            data[current_key] = buffer_lines

    # -------------------------------------------------------------
    # Slicing Array Subsections: [Values, Errors]
    # -------------------------------------------------------------
    # Matrix (Alternating columns: val1, err1, val2, err2, val3, err3)
    matrix_raw = np.array(data["matrix"])
    matrix_vals = matrix_raw[:, 0::2]
    matrix_errs = matrix_raw[:, 1::2]
    
    # Vectors (2 columns: val, err)
    Riso_raw = np.array(data["Riso"])
    Riso_vals = Riso_raw[:, 0]
    Riso_errs = Riso_raw[:, 1]
    
    Rsim_raw = np.array(data["Rsim"])
    Rsim_vals = Rsim_raw[:, 0]
    Rsim_errs = Rsim_raw[:, 1]
    
    sol_raw = np.array(data["solution"])
    sol_vals = sol_raw[:, 0]
    sol_errs = sol_raw[:, 1]

    return (matrix_vals, matrix_errs, 
            Riso_vals, Riso_errs, 
            Rsim_vals, Rsim_errs, 
            sol_vals, sol_errs)

fig, ax = plt.subplots(figsize=(8, 7))
fig1, ax1 = plt.subplots(figsize=(8, 7))
fig2, ax2 = plt.subplots(figsize=(8, 7))
for i, e in enumerate(["B64", "C80", "D96", "E112"]):
    (matrix, matrix_errs, 
        Riso, Riso_errs, 
        Rsim, Rsim_errs, 
        sol, sol_errs) = parse_jackknife_with_errors(f"system_{e}_WP25.txt")
    print(Riso, Rsim, sol)
    # Derive vector steps from column values
    v1 = matrix[:, 0] * sol[0]
    v2 = matrix[:, 1] * sol[1]
    v3 = matrix[:, 2] * sol[2]

    # Chain tracking coordinates
    step1 = Rsim + v1
    step2 = step1 + v2
    LHS = step2 + v3


    # Scatter milestones
    #print(r"$R_{\rm sim}$" + e)
    lab = r"$P_{{\rm sim}} \text{{{}}}$".format(e)
    labL = r"$P^{{\rm WP25'}} \text{{{}}}$".format(e)
      # 2. Assign the label only on the first iteration (i == 0)
    label_iso = r"$P^{{\rm WP25}}$" if e == "B64" else r"_$P^{{\rm WP25}}$"
    print(lab)
    ax.scatter(Riso[0], Riso[1], marker='s', color='navy', s=80, label=label_iso, zorder=5)
    ax.scatter(Rsim[0], Rsim[1], marker='^', s=100, label=lab, zorder=4)
    ax.scatter(LHS[0], LHS[1], marker='o', color='magenta', s=100, label=labL, zorder=4)

    ax.scatter(step1[0], step1[1], marker='s', color='navy', s=80, label='', zorder=5,visible=False)
    ax.scatter(step2[0], step2[1], marker='s', color='navy', s=80, label='', zorder=5,visible=False)
    ax.scatter(LHS[0], LHS[1], marker='s', color='navy', s=80, label='', zorder=5,visible=False)
    # Sequential arrow mapping
    ax.annotate('', xy=step1[0:2], xytext=Rsim[0:2], arrowprops=dict(arrowstyle="->", color="darkorange", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0))
    # ax.plot([], [], color='darkorange', lw=2, label=r'Vector 1: $(M_{11}x_1, M_{21}x_1)$')

    ax.annotate('', xy=step2[0:2], xytext=step1[0:2], arrowprops=dict(arrowstyle="->", color="dodgerblue", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0,linestyle="--"))
    # ax.plot([], [], color='dodgerblue', lw=2, label=r'Vector 2: $(M_{12}x_2, M_{22}x_2)$')

    ax.annotate('', xy=LHS[0:2], xytext=step2[0:2], arrowprops=dict(arrowstyle="->", color="purple", lw=2, mutation_scale=12,shrinkA=0, linestyle=":"))
    # ax.plot([], [], color='purple', lw=2, label=r'Vector 3: $(M_{13}x_3, M_{23}x_3)$')

    # Adjust layout window geometry 
    # plt.tight_layout()
    # plt.tight_layout()
    # plt.savefig("system_B64_WP25_ls.pdf", format="pdf")
    # plt.close()
    ###############################################################################àà 
    ax1.scatter(Riso[1], Riso[2], marker='s', color='navy', s=80, label=label_iso, zorder=5)
    ax1.scatter(Rsim[1], Rsim[2], marker='^', s=100, label=lab, zorder=4)
    ax1.scatter(LHS[1], LHS[2], marker='o', color='magenta', s=100, label=labL, zorder=4)

    ax1.scatter(step1[1], step1[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    ax1.scatter(step2[1], step2[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    ax1.scatter(LHS[1], LHS[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    # Sequential arrow mapping
    ax1.annotate('', xy=step1[1:3], xytext=Rsim[1:3], arrowprops=dict(arrowstyle="->", color="darkorange", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0))
    # ax1.plot([], [], color='darkorange', lw=2, label=r'Vector 1: $(M_{11}x_1, M_{21}x_1)$')

    ax1.annotate('', xy=step2[1:3], xytext=step1[1:3], arrowprops=dict(arrowstyle="->", color="dodgerblue", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0,linestyle="--"))
    # ax1.plot([], [], color='dodgerblue', lw=2, label=r'Vector 2: $(M_{12}x_2, M_{22}x_2)$')

    ax1.annotate('', xy=LHS[1:3], xytext=step2[1:3], arrowprops=dict(arrowstyle="->", color="purple", lw=2, mutation_scale=12,shrinkA=0, linestyle=":"))
    # ax1.plot([], [], color='purple', lw=2, label=r'Vector 3: $(M_{13}x_3, M_{23}x_3)$')

    ###############################################################################àà 
    ax2.scatter(Riso[0], Riso[2], marker='s', color='navy', s=80, label=label_iso, zorder=5)
    ax2.scatter(Rsim[0], Rsim[2], marker='^', s=100, label=lab, zorder=4)
    ax2.scatter(LHS[0], LHS[2], marker='o', color='magenta', s=100, label=labL, zorder=4)

    ax2.scatter(step1[0], step1[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    ax2.scatter(step2[0], step2[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    ax2.scatter(LHS[0], LHS[2], marker='s', color='navy', s=8,  zorder=5,visible=False)
    # Sequential arrow mapping
    ax2.annotate('', xy=step1[[0,2]], xytext=Rsim[[0,2]], arrowprops=dict(arrowstyle="->", color="darkorange", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0))
    # ax1.plot([], [], color='darkorange', lw=2, label=r'Vector 1: $(M_{11}x_1, M_{21}x_1)$')

    ax2.annotate('', xy=step2[[0,2]], xytext=step1[[0,2]], arrowprops=dict(arrowstyle="->", color="dodgerblue", lw=2, mutation_scale=12,shrinkA=0,  shrinkB=0,linestyle="--"))
    # ax1.plot([], [], color='dodgerblue', lw=2, label=r'Vector 2: $(M_{12}x_2, M_{22}x_2)$')

    ax2.annotate('', xy=LHS[[0,2]], xytext=step2[[0,2]], arrowprops=dict(arrowstyle="->", color="purple", lw=2, mutation_scale=12,shrinkA=0, linestyle=":"))
    # ax1.plot([], [], color='purple', lw=2, label=r'Vector 3: $(M_{13}x_3, M_{23}x_3)$')



ax.set_xlabel('Component 1 ($x_1$)')
ax.set_ylabel('Component 2 ($x_2$)')
ax.set_title('Column Basis Vector Decomposition Trajectory')
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left')
fig.savefig("system_WP25_ls.pdf", bbox_inches='tight')

# Adjust layout window geometry 
ax1.set_xlabel('Component 2 ($x_2$)')
ax1.set_ylabel('Component 3 ($x_3$)')
ax1.set_title('Column Basis Vector Decomposition Trajectory')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='upper left')
fig1.savefig("system_WP25_sc.pdf", bbox_inches='tight')


ax2.set_xlabel('Component 1 ($x_1$)')
ax2.set_ylabel('Component 3 ($x_3$)')
ax2.set_title('Column Basis Vector Decomposition Trajectory')
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper right')
fig2.savefig("system_WP25_lc.pdf", bbox_inches='tight')