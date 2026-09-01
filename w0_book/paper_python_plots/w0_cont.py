import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
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
green= "#2CA02C"
colors = [orange,blue,maroon,navy,yellow,lavender]
colors_dis = {"tm":red,"OS":blue}
symbol_dis = {"tm":"^","OS":"v"}
########################################################################################
########################################################################################

def plot_fit(ax, basename, var, data_type=None, noribbon=False,
             id_x=1, noline=False, labelfit="fit", width=0.02, size=1,
             id_color=None, id_shape=None, single_name_for_fit=None,
             nolabel_for_fit=False, nudge=0, alpha_line=1, alpha_ribbon=0.5,
             stroke=1, filter_data=None):
    """
    Translates the R plot_fit logic into Matplotlib.
    Instead of passing a 'gg' object, we pass a Matplotlib Axis object ('ax').
    """
    filed = f"{basename}_fit_data.txt"
    if not os.path.exists(filed):
        raise FileNotFoundError(f"CRITICAL ERROR: The file '{filed}' does not exist in the current directory: {os.getcwd()}")
        
    # Read text data (handling white-space separation)
    df = pd.read_csv(filed, header=None, sep=r'\s+')
    
    if filter_data is not None:
        last_col_idx = df.shape[1] - 1
        df = df[df.iloc[:, last_col_idx].isin(filter_data)].copy()
        
    # Python 0-based index vs R 1-based index translation
    idy_col = df.shape[1] - 3  # R logic: ncol - 2 becomes python index - 3
    idx_col = id_x - 1
    
    # Track the final grouping tracking columns
    fit_col_idx = df.shape[1] - 1
    
    color_col_idx = fit_col_idx if id_color is None else (id_color - 1)
    shape_col_idx = fit_col_idx if id_shape is None else (id_shape - 1)
    
    color_type = df.iloc[:, color_col_idx].astype(str).tolist()
    shape_type = df.iloc[:, shape_col_idx].astype(str).tolist()
    
    start_fit = int(df.iloc[0, fit_col_idx])
    end_fit = int(df.iloc[-1, fit_col_idx])
    Nfits = list(range(start_fit, end_fit + 1))
    
    if data_type is not None:
        if isinstance(data_type, (str, int, float)) or len(data_type) == 1:
            val = data_type[0] if isinstance(data_type, list) else data_type
            color_type = [str(val)] * len(df)
            shape_type = [str(val)] * len(df)
        else:
            color_type = df.iloc[:, fit_col_idx].tolist()
            shape_type = df.iloc[:, fit_col_idx].tolist()
            count = 0
            for n in range(len(df)):
                if n != 0 and df.iloc[n, fit_col_idx] != df.iloc[n - 1, fit_col_idx]:
                    count += 1
                color_type[n] = str(data_type[count])
                shape_type[n] = str(data_type[count])
                
    mycol = [f"{labelfit}{c}" for c in set(color_type)]
    if data_type is not None:
        mycol = list(Nfits)
        count = 0
        mycol[0] = f"{labelfit}{data_type[0]}"
        for n in range(len(df)):
            if n != 0 and df.iloc[n, fit_col_idx] != df.iloc[n - 1, fit_col_idx]:
                count += 1
                mycol[count] = f"{labelfit}{data_type[count]}"
                
    if len(mycol) != len(Nfits):
        mycol = [f"{labelfit}{n}" for n in Nfits]
        
    if single_name_for_fit is not None:
        mycol = [single_name_for_fit] * len(Nfits)

    # 1. Plot the fits (Ribbons and Lines)
    if (not noribbon) or (not noline):
        for idx, n in enumerate(Nfits):
            file = f"{basename}_fit_out_n{n}_{var}.txt"
            if not os.path.exists(file):
                raise FileNotFoundError(f"CRITICAL ERROR: The file '{file}' does not exist.")
            
            sub_df = pd.read_csv(file, header=None, sep=r'\s+')
            sub_df.columns = ['x', 'fit', 'fiterr']
            
            x_vals = sub_df['x'] + nudge
            y_vals = sub_df['fit']
            ymin = y_vals - sub_df['fiterr']
            ymax = y_vals + sub_df['fiterr']
            lbl = str(mycol[idx])
            
            if not noribbon:
                ax.fill_between(x_vals, ymin, ymax, alpha=alpha_ribbon, label=lbl)
            if not noline:
                # ax.plot(x_vals, y_vals, alpha=alpha_line, label=lbl, color="gray", linewidth=0.5)
                # ax.plot(x_vals, y_vals,color="black",alpha=0.2,linewidth=0.5)
                ax.plot(x_vals, y_vals,color=green, alpha=0.5,linestyle="dashed",linewidth=0.4)

    # 2. Plot Points & Errorbars from main text data frame
    # Collect data arrays
    x_data = df.iloc[:, idx_col] + nudge
    y_data = df.iloc[:, idy_col]
    y_err = df.iloc[:, idy_col + 1]
    
    # We dynamically handle categorical colors/shapes by splitting scatter plots
    unique_combos = sorted(list(set(zip(color_type, shape_type))))
    
    # Marker map variant safely addressing varied markers
    # marker_choices = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*']
    marker_choices = ['d', 'v', '^', 'D', 's', '<', '>', 'p', '*']
    marker_map = {combo: marker_choices[i % len(marker_choices)] for i, combo in enumerate(unique_combos)}
    #marker_map = ["o","^","^"]

    for combo in unique_combos:
        mask = [c == combo[0] and s == combo[1] for c, s in zip(color_type, shape_type)]
        if any(mask):
            ax.errorbar(
                x_data[mask], y_data[mask], yerr=y_err[mask],
                fmt=marker_map[combo],
                color=green,
                elinewidth=size, capsize=width*1000, 
                # markersize=size*5,
                # markeredgewidth=stroke,
                label=None
            )
def read_fit_file(file_path):
    # Read space-separated data (equivalent to read.table with fill=True)
    # Generate 40 columns (0 to 39 in Python's 0-based indexing)
    df = pd.read_csv(file_path, sep=r'\s+', header=None, names=range(40))
    
    # Store results in a Python dictionary (equivalent to R's list)
    result = {}
    
    # R: df[1,2] -> Row 1, Col 2 (1-based index) -> Python: .iloc[0, 1]
    result['npar'] = int(df.iloc[0, 1])
    result['chi2dof'] = float(df.iloc[1, 1])
    
    npar = result['npar']
    
    # R: Row range c(3 : 3+npar-1) -> Python slice: 2 to 2+npar
    # R: Column range c(1,2,3) -> Python column indices: 0, 1, 2
    result['P'] = df.iloc[2 : 2 + npar, [0, 1, 2]].reset_index(drop=True)
    
    # R: Row range c(3+npar : 3+npar*2-1) -> Python slice: 2+npar to 2+npar*2
    # R: Column range c(1:npar) -> Python column range: 0 to npar
    result['C'] = df.iloc[2 + npar : 2 + npar * 2, list(range(npar))].reset_index(drop=True)
    
    # R: df[3+npar*2, 2] -> Python offset index: 2 + npar*2
    result['dof'] = int(df.iloc[2 + npar * 2, 1])
    
    # R: df[3+npar*2+1, 2] -> Python offset index: 3 + npar*2
    result['ndata'] = int(df.iloc[3 + npar * 2, 1])
    
    return result

def calculate_baic_average(v, err, chi2dof, dof, npar, multiplicity=1):
    # Convert inputs to numpy arrays to ensure vector math works correctly
    v = np.asarray(v)
    err = np.asarray(err)
    chi2dof = np.asarray(chi2dof)
    dof = np.asarray(dof)
    npar = np.asarray(npar)
    
    result = {}
    
    # Calculate degrees of freedom + parameters
    Nmeas = dof + npar
    
    # R: AIC <- exp(-0.5 * (chi2dof * dof + 2 * npar - 2 * Nmeas)) / multiplicity
    aic = np.exp(-0.5 * (chi2dof * dof + 2 * npar - 2 * Nmeas)) / multiplicity
    
    # Normalize AIC weights
    N = np.sum(aic)
    aic = aic / N
    result['AIC'] = aic
    
    # Calculate the weighted mean
    result['m'] = np.sum(v * aic)
    
    # Calculate statistical and systematic variance components
    stat_variance = np.sum(aic * (err ** 2))
    syst_variance = np.sum(aic * ((v - result['m']) ** 2))
    
    # Store standard deviations (square roots of variances)
    result['dm'] = np.sqrt(stat_variance + syst_variance)
    result['stat'] = np.sqrt(stat_variance)
    result['syst'] = np.sqrt(syst_variance)
    
    return result

# --- Main Script Execution ---

C = -5
path = "/home/garofalo/analysis/flow/data_20/fit_all_beta/"
basenames = [
    f"fit_w0_lin_deriv_FLAG_a2",
    f"fit_w0_lin_deriv_FLAG_a2_noB",
    f"fit_w0_lin_deriv_FLAG_a2_noC",
    f"fit_w0_lin_deriv_FLAG_a2_noD",
    f"fit_w0_lin_deriv_FLAG_a2_noE",
    f"fit_w0_lin_deriv_FLAG_a2_a4",
    # f"fit_w0_lin_deriv_FLAG_a2_Husung0.42"
]

count = len(basenames)

df = pd.DataFrame({
    "fit": [""] * count,
    "res": [0.0] * count,
    "err": [0.0] * count,
    "chi2dof": [0.0] * count,
    "dof": [0] * count,
    "Npar": [0] * count,
    "Ndat": [0] * count,
    "mult": [0.0] * count
    })

# Generate labels, cleaning syntax patterns natively
legend_name = [re.sub(r"fit_fpi_|\.000000", "", name) for name in basenames]
# legend_name = [f"\\verb|{name}|" for name in legend_name]

labels = [["WTI", "tm", "OS"]]

# Initialize the Matplotlib figure canvas
# Defaulting layout variables width/height if missing in original snippet scope
width, height = 800, 600 
#fig, ax = plt.subplots(figsize=(width / 100, height / 100))
fig, ax = plt.subplots(1, 2, width_ratios=[1.,2.], sharey=True,figsize=(12,6))

Nboot=2000
flat_boot = [] 
# Iterate and append layers directly onto the initialized axes
for j, basename in enumerate(basenames):
    plot_fit(
        ax=ax[1],
        basename=os.path.join(path, basename),
        var="a2",
        data_type=labels[0],
        id_x=1,
        single_name_for_fit="",
        width=0.004,
        size=0.8,
        nudge=0,
        noline=False,
        noribbon=True,
        alpha_line = 0.5,
        stroke=0.1
    )
    # Assuming 'path', 'basenames', and 'j' are defined in your loop:
    file_path = os.path.join(path, f"{basenames[j]}_fit_P.dat")

    # Call the converted python function
    fit = read_fit_file(file_path)
    # Initialize row j with a list of values
    # Note: The list must have the exact same number of elements as there are columns (9 columns)
    df.iloc[j] = [basenames[j], fit['P'].iloc[0,1], fit['P'].iloc[0,2], fit['chi2dof'], fit['dof'], fit['npar'], fit['ndata'], 1]
    bt=np.random.normal(fit['P'].iloc[0,1], fit['P'].iloc[0,2], Nboot)
    flat_boot.extend(bt)
    
ave_BAIC = calculate_baic_average(
    v=df['res'].to_numpy(),
    err=df['err'].to_numpy(),
    chi2dof=df['chi2dof'].to_numpy(),
    dof=df['dof'].to_numpy(),
    npar=df['Npar'].to_numpy(),
    multiplicity=df['mult'].to_numpy()
)
BAIC_ave = ave_BAIC['m']
BAIC_err = ave_BAIC['dm']
stat_err = ave_BAIC['stat']
flat_weig =[]
for j, basename in enumerate(basenames):
    flat_weig.extend([ave_BAIC['AIC'][j]]*Nboot)

# Reference benchmark flag lines
fpi_FLAG = 0.17236
ax[1].axhline(y=fpi_FLAG, color='black', linestyle='--', label='wp25')
ax[1].scatter([0], [fpi_FLAG], color='red', marker='x', label='wp25')

# Typography and Axis setup
title = ""
xlabel = r"$a^2$ [fm$^2$]"
ylabel = r"$w_{0}$ [fm]"

legend_position = (0.7, 0.98)

if title != "":
    ax[1].set_title(title)
if xlabel != "":
    ax[1].set_xlabel(xlabel)
if ylabel != "":
    ax[0].set_ylabel(ylabel)

# Customizing the Classic Matplotlib border presentation (Replicates theme_matplotlib)
ax[1].patch.set_facecolor('white')
for spine in ax[1].spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)
# ax[1].tick_params(colors='black', direction='out')
ax[1].yaxis.get_label().set_visible(False)
ax[1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
# Clean duplicate handles generated during the iterative subplot loops
handles, plot_labels = ax[1].get_legend_handles_labels()
by_label = dict(zip(plot_labels, handles))
# ax[1].legend(by_label.values(), by_label.keys(), loc="upper left",
#             #  bbox_to_anchor=legend_position
#              )
legend(ax[1],by_label.values(), by_label.keys(), loc="upper left",
            #  bbox_to_anchor=legend_position
             )

ax[1].errorbar([0.],[BAIC_ave],[BAIC_err],fmt="x",color="black",label=f"BAIC average")
ax[1].errorbar([0.],[BAIC_ave],[stat_err],fmt="",color="black")
ax[1].set_ylim(0.167, 0.176)
# histogram
# Nota per Marco: in flat_boot c'è sostanzialmente un array che contiene un ricampionamento bootstrap
# dei limiti al continuo. flat_weig sono i pesi normalizzati con la somma dei pesi (e li moltiplico
# per 100 per ottenere la percentuale). Ogni modello ha il suo peso che è lo stesso per ogni bootstrap.
# Nota importante: flat_boot ha dimensione (numero di modelli)*(numero di samples) ma è un array 1D. Lo stesso
# vale per flat_weights, dove ripeti (numero di samples) volte lo stesso weight per ogni modello.
ax[0].hist(flat_boot, bins=20, weights=np.array(flat_weig)*100/Nboot,orientation="horizontal",histtype='step', linewidth=1.5,color="black")
print(np.array(flat_weig).sum()/Nboot )
ax[0].set_xlim(ax[0].get_xlim())
# BAIC_ave e BAIC_err si capisce cosa sono e te li devi calcolare a parte
ax[0].fill_between(ax[0].get_xlim(),[BAIC_ave-BAIC_err,BAIC_ave-BAIC_err],
                            [BAIC_ave+BAIC_err,BAIC_ave+BAIC_err],color=lavender,alpha=0.2)
ax[0].fill_between(ax[0].get_xlim(),[BAIC_ave-stat_err,BAIC_ave-stat_err],
                            [BAIC_ave+stat_err,BAIC_ave+stat_err],color=lavender,alpha=0.2)

ax[0].set_xlabel(r"%")
ax[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax[0].grid(True, which='minor', axis='x')
ax[0].tick_params(axis='x',which='minor',size=0)

plt.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.12, wspace=0)

# Save configuration
# Matplotlib saves vector figures cleanly via .pdf or .svg. 
# If your final step compiles in LaTeX via pgf/tikz, use .pgf extension format target instead.


fpi3reg = "w0cont"
# plt.tight_layout()
plt.savefig(f"{fpi3reg}.pdf", format="pdf")
plt.close()
