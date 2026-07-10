import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import FormatStrFormatter
import re
import numpy as np

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



DEFAULT_PROJECTIONS = [
    ("R_pi", "R_Ds"),
    ("R_K", "R_Ds"),
    ("R_pi", "R_K"),
]
DEFAULT_LINEAR_ORDER = [1, 2, 0]
DEFAULT_QUARKS = ["l", "s", "c"]
DEFAULT_COLORS = {"l": "blue", "s": "green", "c": "red"}
DEFAULT_LINESTYLES = {"l": "-", "s": "--", "c": ":"}
DEFAULT_LABELS = {
    "R_pi": r"$R_\pi$",
    "R_K": r"$R_K$",
    "R_Ds": r"$R_{D_s}$",
}


def configure_mistuning_plot_style():
    plt.rcParams.update({
        "text.usetex": False,
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "font.family": "serif",
        "font.serif": ["cmr10"],
        "font.size": 18,
        "axes.labelsize": 18,
        "axes.titlesize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "figure.titlesize": 18,
        "errorbar.capsize": 5,
        "lines.markeredgewidth": 2.0,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "axes.grid": True,
        "axes.grid.which": "major",
        "grid.color": "gray",
        "grid.linestyle": "--",
        "grid.linewidth": 0.5,
        "grid.alpha": 0.7,
    })


def draw_styled_arrow(ax, x, y, dx, dy, color, linestyle):
    head_fraction = 0.05
    shaft_end_x = x + (1.0 - head_fraction) * dx
    shaft_end_y = y + (1.0 - head_fraction) * dy
    ax.add_line(Line2D(
        [x, shaft_end_x],
        [y, shaft_end_y],
        linewidth=2.4,
        linestyle=linestyle,
        color=color,
        transform=ax.transData,
        zorder=9000,
    ))
    arrow_head = FancyArrowPatch(
        (shaft_end_x, shaft_end_y),
        (x + dx, y + dy),
        arrowstyle="simple",
        mutation_scale=16,
        linewidth=2.4,
        linestyle="-",
        color=color,
        shrinkA=0,
        shrinkB=0,
        transform=ax.transData,
        zorder=9000,
    )
    ax.add_patch(arrow_head)
    return arrow_head


def draw_arrow_key(ax, x, y, length, label, color, linestyle, coordinates="axes"):
    if coordinates == "axes":
        transform = ax.transAxes
    elif coordinates == "figure":
        transform = ax.figure.transFigure
    else:
        transform = ax.transData
    ax.add_line(Line2D(
        [x, x + length],
        [y, y],
        linewidth=2.4,
        linestyle=linestyle,
        color=color,
        transform=transform,
        clip_on=False,
        zorder=10000,
    ))
    key_arrow_head = Line2D(
        [x + length],
        [y],
        marker=">",
        markersize=8,
        linestyle="",
        markerfacecolor=color,
        markeredgecolor=color,
        color=color,
        transform=transform,
        clip_on=False,
        zorder=10000,
    )
    ax.add_line(key_arrow_head)
    ax.text(
        x + length + 0.014,
        y,
        label,
        transform=transform,
        va="center",
        ha="left",
        fontsize=plt.rcParams["legend.fontsize"],
        clip_on=False,
        zorder=10000,
    )
    return key_arrow_head


def draw_marker_key(ax, x, y, marker, label, color, coordinates="figure"):
    transform = ax.figure.transFigure if coordinates == "figure" else ax.transAxes
    ax.add_line(Line2D(
        [x],
        [y],
        marker=marker,
        markersize=14,
        linestyle="",
        markerfacecolor=color,
        markeredgecolor=color,
        transform=transform,
        clip_on=False,
        zorder=10000,
    ))
    ax.text(
        x + 0.014,
        y,
        label,
        transform=transform,
        va="center",
        ha="left",
        fontsize=plt.rcParams["legend.fontsize"],
        clip_on=False,
        zorder=10000,
    )


def _point_xy(point, x_obs, y_obs):
    return point[x_obs], point[y_obs]


def _vector_xy(vectors, projection_index, x_obs, y_obs):
    vector = vectors[projection_index]
    if isinstance(vector, dict):
        return vector[x_obs], vector[y_obs]
    return vector[0], vector[1]
# Map string keys to numerical indices for sequence handling fallback
# OBS_INDEX_MAP = {"R_pi": 0, "R_K": 1, "R_Ds": 2}

# def _point_xy(point, x_obs, y_obs):
#     """
#     Safely extracts coordinates whether point is a dict or an indexed sequence.
#     """
#     if isinstance(point, dict):
#         return point[x_obs], point[y_obs]
    
#     # Fallback for arrays/lists using the mapping
#     idx_x = OBS_INDEX_MAP[x_obs]
#     idx_y = OBS_INDEX_MAP[y_obs]
#     return point[idx_x], point[idx_y]


# def _vector_xy(vectors, projection_index, x_obs, y_obs):
#     """
#     Safely extracts vector deltas from dicts, nested lists, or flat slices.
#     """
#     vector = vectors[projection_index]
#     if isinstance(vector, dict):
#         return vector[x_obs], vector[y_obs]
        
#     # If it is a nested collection of vectors per projection
#     if hasattr(vector, '__len__') and len(vector) >= 2:
#         # If it's a flat 2-element vector coordinate pair (dx, dy)
#         if not hasattr(vector[0], '__len__'):
#             return vector[0], vector[1]
#         # If it's tracking specific mapped spatial slices
#         idx_x = OBS_INDEX_MAP[x_obs]
#         idx_y = OBS_INDEX_MAP[y_obs]
#         return vector[idx_x], vector[idx_y]
        
#     return vector, vector

def _ensemble_label(ensemble):
    return ensemble[0] if isinstance(ensemble, str) else str(ensemble)


def plot_mistuning_planes(
    ensembles,
    origins,
    vectors,
    flag_points,
    ax_lim,
    iso_point,
    labels=None,
    projections=None,
    linear_order=None,
    save_ortogonal_path=None,
):
    """Draw the linear and orthogonal-projection mistuning plots.

    Parameters
    ----------
    ensembles:
        Iterable of ensemble names. The first character is used as the plotted label.
    origins:
        Dict mapping ensemble name to observable coordinates, e.g.
        {"B64": {"R_pi": ..., "R_K": ..., "R_Ds": ...}, ...}.
    vectors:
        Dict mapping ensemble name -> quark -> list of three projected vectors.
        Each projected vector can be either ``(dx, dy)`` or a dict with the
        projection observable names as keys. Example:
        {"B64": {"l": [(dx0, dy0), (dx1, dy1), (dx2, dy2)], ...}, ...}.
    flag_points:
        Observable coordinates of the FLAG-scheme point, or a list of three
        per-projection points.
    ax_lim:
        Dict mapping observable name to ``[min, max]`` limits.
    labels:
        Optional dict mapping observable names to axis labels.
    projections:
        Optional list of three ``(x_obs, y_obs)`` pairs.
    linear_order:
        Optional list mapping projection index to linear-panel index.
    save_ortogonal_path:
        If not None, save the orthogonal-projection figure to this path.

    Returns
    -------
    fig, ax, ffigg, axs
        The linear figure/axes and orthogonal-projection figure/axes.
    """
    configure_mistuning_plot_style()
    labels = DEFAULT_LABELS if labels is None else labels
    projections = DEFAULT_PROJECTIONS if projections is None else projections
    linear_order = DEFAULT_LINEAR_ORDER if linear_order is None else linear_order

    ffigg = plt.figure(figsize=(8.4, 8))
    gs = ffigg.add_gridspec(
        2,
        2,
        width_ratios=[1, 1],
        height_ratios=[1, 1],
        top=0.98,
        bottom=0.08,
        left=0.12,
        right=0.98,
        hspace=.05,
        wspace=.05,
    )
    ax1 = ffigg.add_subplot(gs[0, 0])
    ax2 = ffigg.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = ffigg.add_subplot(gs[1, 0], sharex=ax1)
    axs = [[ax1, ax2], [ax3]]

    fig, ax = plt.subplots(1, 3, figsize=(14, 5.4))
    fig.subplots_adjust(left=0.06, right=0.995, bottom=0.12, top=0.85, wspace=0.04)

    for pp, (x_obs, y_obs) in enumerate(projections):
        linear_ax = ax[linear_order[pp]]
        ortho_ax = axs[pp // 2][pp % 2]
        linear_ax.set_box_aspect(4 / 3.5)
        ortho_ax.set_box_aspect(1)

        for idx, ens in enumerate(ensembles):
            _draw_ensemble_label(linear_ax, origins[ens], ens, x_obs, y_obs, ax_lim, pp, idx, linear=True)
            _draw_ensemble_label(ortho_ax, origins[ens], ens, x_obs, y_obs, ax_lim, pp, idx, linear=False)

            for q in DEFAULT_QUARKS:
                dx, dy = _vector_xy(vectors[ens][q], pp, x_obs, y_obs)
                draw_styled_arrow(
                    linear_ax,
                    origins[ens][x_obs],
                    origins[ens][y_obs],
                    dx,
                    dy,
                    DEFAULT_COLORS[q],
                    DEFAULT_LINESTYLES[q],
                )
                draw_styled_arrow(
                    ortho_ax,
                    origins[ens][x_obs],
                    origins[ens][y_obs],
                    dx,
                    dy,
                    DEFAULT_COLORS[q],
                    DEFAULT_LINESTYLES[q],
                )
                if idx == 0 and pp == 2:
                    quark_fancy = r"\ell" if q == "l" else q
                    draw_arrow_key(
                        ortho_ax,
                        1.42,
                        0.585 - DEFAULT_QUARKS.index(q) * 0.15,
                        0.13,
                        rf"$\vec{{r}}_{quark_fancy}$",
                        DEFAULT_COLORS[q],
                        DEFAULT_LINESTYLES[q],
                    )

        _format_projection_axis(linear_ax, x_obs, y_obs, labels, ax_lim)
        _format_projection_axis(ortho_ax, x_obs, y_obs, labels, ax_lim)

        my_marker = ['o','s','D','X']
        
        for i, e in enumerate(ensembles):
            print(flag_points[e])
            flag_point = flag_points[e][pp] if isinstance(flag_points[e], (list, tuple)) else flag_points[e]
            flag_x, flag_y = _point_xy(flag_point, x_obs, y_obs)
            linear_ax.scatter([flag_x], [flag_y],
                            marker=my_marker[i],
                              s=100, 
                              color="black"
                              )
            ortho_label = f"WP25' scheme {e}" if pp == 2 else ""
            ortho_ax.scatter([flag_x], [flag_y],
                            #  marker="*",
                            #  s=200, 
                             color="black",
                             label=ortho_label)
        flag_point = iso_point[pp] if isinstance(iso_point, (list, tuple)) else iso_point
        flag_x, flag_y = _point_xy(iso_point, x_obs, y_obs)
        linear_ax.scatter([flag_x], [flag_y], marker="*", s=200, color="black")
        ortho_label = f"WP25 scheme" if pp == 2 else ""
        ortho_ax.scatter([flag_x], [flag_y], marker="*", s=200, 
                        color="black",
                            label=ortho_label)
    _draw_linear_legend(fig, ax[0],my_marker,ensembles)
    _draw_orthogonal_legend(axs, ax_lim)

    if save_ortogonal_path is not None:
        ffigg.savefig(save_ortogonal_path)

    return fig, ax, ffigg, axs


def _draw_ensemble_label(ax, point, ensemble, x_obs, y_obs, ax_lim, projection_index, ensemble_index, linear):
    if linear:
        directions = [
            [(2, 2), (2, 5.4), (2, 2), (2, -5.4)],
            [(-1, -6), (0, -6), (-1, -6), (3, -5.)],
            [(3, 2), (3, 2), (3, 2), (-5.5, -5.5)],
        ]
    else:
        directions = [
            [(1.5, -7), (3, 7), (3, -6), (4, -7)],
            [(-1, -8), (0, -8), (-3, -7), (-3, -7.)],
            [(4, 1), (3, 2), (4, 1), (-8, -8)],
        ]
    direction = directions[projection_index][ensemble_index % len(directions[projection_index])]
    dx = direction[0] * 0.01 * (ax_lim[x_obs][1] - ax_lim[x_obs][0])
    dy = direction[1] * 0.01 * (ax_lim[y_obs][1] - ax_lim[y_obs][0])
    ax.text(
        point[x_obs] + dx,
        point[y_obs] + dy,
        _ensemble_label(ensemble),
        bbox=dict(facecolor="white", edgecolor="gray", alpha=0.6, boxstyle="round,pad=0.2"),
        fontsize=18,
        zorder=10000,
    )


def _format_projection_axis(ax, x_obs, y_obs, labels, ax_lim):
    ax.set_xlabel(labels[x_obs])
    ax.set_xlim(ax_lim[x_obs])
    ax.set_ylabel(labels[y_obs])
    ax.set_ylim(ax_lim[y_obs])


def _draw_linear_legend(fig, anchor_ax,my_marker,ensembles):
    y_box = 0.88
    height_labels = 0.92
    fig.patches.append(FancyBboxPatch(
        (0.1, y_box),
        0.85,
        0.09,
        boxstyle="round,pad=0.012",
        transform=fig.transFigure,
        facecolor="none",
        edgecolor="gray",
        linewidth=1.0,
        zorder=9990,
    ))
    for id_q, q in enumerate(DEFAULT_QUARKS):
        quark_fancy = r"\ell" if q == "l" else q
        draw_arrow_key(
            anchor_ax,
            0.1 + id_q * 0.1,
            height_labels,
            0.045,
            rf"$\vec{{r}}_{quark_fancy}$",
            DEFAULT_COLORS[q],
            DEFAULT_LINESTYLES[q],
            coordinates="figure",
        )
    start_x = 0.1 + 3 * 0.1
    start_height = height_labels+0.03
    current_height = start_height
    draw_marker_key(anchor_ax, start_x, start_height, "*", "WP25 scheme", "black")
    line_spacing = 0.05  # Adjust this value to make the gap larger or smaller
    current_height = start_height - ( line_spacing)
    for i, e in enumerate(ensembles):
        draw_marker_key(anchor_ax, start_x, current_height, my_marker[i], f"WP25' scheme {e}", "black")
        if (i%2==0):
            current_height = start_height
            start_x += 0.2
        if(i==1):
            current_height = start_height - ( line_spacing)


def _draw_orthogonal_legend(axs, ax_lim):
    handles, labels = axs[1][0].get_legend_handles_labels()
    dummy = Line2D([], [], linestyle="", alpha=0)
    handles.extend([dummy, dummy, dummy])
    labels.extend(["", "", ""])
    leg = axs[1][0].legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.2, 0.84),
        borderpad=1.2,
        labelspacing=1.62,
    )
    leg.set_zorder(1)

    axs[0][1].set_xlim(ax_lim["R_K"])
    axs[0][1].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ticks_K = axs[0][1].get_xticks()
    axs[1][0].set_yticks(ticks_K)
    axs[1][0].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axs[1][0].set_ylim(ax_lim["R_K"])

    axs[0][1].invert_xaxis()
    axs[0][1].set(xlabel=None)
    axs[0][0].set(xlabel=None)
    axs[0][0].tick_params(axis="x", labelbottom=False)
    axs[0][1].set(ylabel=None)
    axs[0][1].tick_params(axis="y", labelleft=False)

ensembles = ["B64", "C80", "D96", "E112"]
origins = {}
vectors = {}
wp25_points = {}
wp251_points = {}
# vector order
# 0: ("R_pi", "R_Ds")
# 1: ("R_K", "R_Ds")
# 2: ("R_pi", "R_K")
# origins = {
#     "B64": {"R_pi": 1.06, "R_K": 14.0, "R_Ds": 14.9},
#     "C80": {"R_pi": 1.08, "R_K": 14.2, "R_Ds": 15.0},
# }
ax_lim = {
    "R_pi": [0.116, 0.123],
    "R_K": [0.415, 0.44],
    "R_Ds": [1.67, 1.72]
}
for i, e in enumerate(ensembles):
    (matrix, matrix_errs, 
        Riso, Riso_errs, 
        Rsim, Rsim_errs, 
        sol, sol_errs) = parse_jackknife_with_errors(f"system_{e}_WP25.txt")
    print(Riso, Rsim, sol)
    origins[e] = {
        "R_pi": Rsim[0],
        "R_K": Rsim[1],
        "R_Ds": Rsim[2]
    }
    v1 = matrix[:, 0] * sol[0]
    v2 = matrix[:, 1] * sol[1]
    v3 = matrix[:, 2] * sol[2]
    vectors[e] = {
        "l" : [(v1[0], v1[2]), (v1[1], v1[2]), (v1[0], v1[1])],
        "s" : [(v2[0], v2[2]), (v2[1], v2[2]), (v2[0], v2[1])],
        "c" : [(v3[0], v3[2]), (v3[1], v3[2]), (v3[0], v3[1])]
    }
    step1 = Rsim + v1
    step2 = step1 + v2
    LHS = step2 + v3
    wp251_points[e] = {
        "R_pi": LHS[0],
        "R_K": LHS[1],
        "R_Ds": LHS[2]
    }
    wp25_points[e] = {
        "R_pi": Riso[0],
        "R_K": Riso[1],
        "R_Ds": Riso[2]
    } 
    
    # print(Riso, Rsim, sol)

fig, ax, ffigg, axs = plot_mistuning_planes(
    ensembles=ensembles,
    origins=origins,
    vectors=vectors,
    flag_points=wp251_points,
    ax_lim=ax_lim,
    iso_point= wp25_points[e]
    
)

fig.savefig("mistuning_wp25_linear.pdf")
ffigg.savefig("mistuning_wp25_orthogonal.pdf")
#plt.show()