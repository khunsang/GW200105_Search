"""Ratio of the sensitive spacetime volumes between the eccentric search and the quasi-circular search,
and fraction of recovered injections in both searches as a function of IFAR threshold.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from cycler import cycler
from pycbc import sensitivity
from pycbc import conversions as conv
from pycbc.io.hdf import HFile
from matplotlib import rcParams

physrev_style = {
    'axes.prop_cycle': cycler('color', ['#5d81b4', '#e09b24', '#8eb031', '#eb6235', '#8678b2', '#c46e1a', '#5c9dc7', '#ffbf00', '#a5609c']) + cycler('ls', ['-', '--', '-.', (0, (1,0.85)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), (0, (5, 1)), ':', (4, (10, 3))]),
    'axes.titlesize': 28,
    'axes.labelsize': 24,
    'axes.formatter.use_mathtext': True,
    'axes.linewidth': 1.75,
    'grid.linewidth': 0.5,
    'grid.linestyle': 'dashed',
    'grid.color': 'xkcd:light gray',
    'lines.linewidth': 1,
    'lines.markersize': 10,
    'hatch.linewidth': 0.25,
    'patch.antialiased': True,
    'patch.force_edgecolor': True,
    'xtick.top': True,
    'xtick.bottom': True,
    'xtick.major.size': 3.0,
    'xtick.minor.size': 1.5,
    'xtick.major.width': 0.5,
    'xtick.minor.width': 0.5,
    'xtick.direction': 'out',
    'xtick.minor.visible': True,
    'xtick.major.top': True,
    'xtick.major.bottom': True,
    'xtick.minor.top': True,
    'xtick.minor.bottom': True,
    'xtick.major.pad': 5.0,
    'xtick.minor.pad': 5.0,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'ytick.left': True,
    'ytick.right': True,
    'ytick.major.size': 3.0,
    'ytick.minor.size': 1.5,
    'ytick.major.width': 0.5,
    'ytick.minor.width': 0.5,
    'ytick.direction': 'out',
    'ytick.minor.visible': True,
    'ytick.major.left': True,
    'ytick.major.right': True,
    'ytick.minor.left': True,
    'ytick.minor.right': True,
    'ytick.major.pad': 5.0,
    'ytick.minor.pad': 5.0,
    'legend.frameon': True,
    'legend.fontsize': 18,
    'legend.handlelength': 1.375,
    'legend.labelspacing': 0.4,
    'legend.columnspacing': 1,
    'legend.facecolor': 'white',
    'legend.edgecolor': 'white',
    'legend.framealpha': 1,
    'legend.title_fontsize': 8,
    'figure.figsize': (4.5, 3.4),
    'figure.subplot.left': 0.125,
    'figure.subplot.bottom': 0.175,
    'figure.subplot.top': 0.95,
    'figure.subplot.right': 0.95,
    'figure.autolayout': False,
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts}',
    'mathtext.fontset': 'stixsans',
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.sans-serif': ['Bitstream Vera Sans'],
    'font.size': 24,
    'path.simplify': True,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'savefig.dpi': 400,
    'backend': 'pdf'
}

plt.rcParams.update(physrev_style)

def log_minor_fmt(x, pos, show_divs=[2, 4, 7]):
    '''
    Function to only show certain minor ticks on the log scale, to avoid overcrowding the y-axis.

    Parameters:
    x: the minor tick value
    pos: the position of the minor tick (required by FuncFormatter, but not used here)
    show_divs: the list of divisors of the decade to show (e.g.,
    '''
    if x == 0:
        return ""
    elif np.isclose(x/10**np.floor(np.log10(x)), show_divs).any():
        return f"{x:g}"
    else:
        return ""


def vratio_error(v_ecc, v_qc, v_ecc_err, v_qc_err):
    """
    Calculate the propagated error for the ratio v_ecc / v_qc using error propagation.

    Parameters:
    v_ecc: The value of the eccentric volume.
    v_qc: The value of the quasi-circular volume.
    v_ecc_err: The error in v_ecc.
    v_qc_err: The error in v_qc.

    Returns:
    float: The propagated error in the ratio.
    """
    vr = v_ecc/v_qc
    vr_error = vr*np.sqrt( (v_qc_err/v_qc)**2 + (v_ecc_err/v_ecc)**2 )
    return vr_error


# Hardcoded variables
eccentric_injection_file = 'H1L1V1_INJ_FIND_ECCENTRIC_INJECTIONS_ECCENTRIC_SEARCH.hdf'
quasicircular_injection_file = 'H1L1V1_INJ_FIND_ECCENTRIC_INJECTIONS_QUASICIRCULAR_SEARCH.hdf'
vt_integration_bin = 200

ifar_thresholds = np.logspace(-1, 4, 100) #in years


emax = 0
emin = .99

with HFile(eccentric_injection_file, 'r') as f:
    # Get the found (at any FAR)/missed injection indices
    eccentric_search_foundi = f['found_after_vetoes/injection_index'][:]
    eccentric_search_missedi = f['missed/after_vetoes'][:]
    # retrieve injection parameters
    eccentric_search_dist = f['injections/distance'][:]
    eccentric_search_exclusive_ifar = f['found_after_vetoes/ifar_exc'][:]
    eccentric_search_e0 = f['injections/eccentricity'][:]
    found_eccentric_distance = eccentric_search_dist[eccentric_search_foundi]
    missed_eccentric_distance = eccentric_search_dist[eccentric_search_missedi]

    emax = max(emax, max(eccentric_search_e0))
    emin = min(emin, min(eccentric_search_e0))


with HFile(quasicircular_injection_file, 'r') as f:
    # Get the found (at any FAR)/missed injection indices
    quasicircular_search_foundi = f['found_after_vetoes/injection_index'][:]
    quasicircular_search_missedi = f['missed/after_vetoes'][:]
    # retrieve injection parameters
    quasicircular_search_dist = f['injections/distance'][:]
    quasicircular_search_exclusive_ifar = f['found_after_vetoes/ifar_exc'][:]
    quasicircular_search_e0 = f['injections/eccentricity'][:]

    #print ("quasicircular_search_exclusive_ifar", quasicircular_search_exclusive_ifar)
    found_quasicircular_e0 = quasicircular_search_e0[quasicircular_search_foundi]
    found_quasicircular_distance = quasicircular_search_dist[quasicircular_search_foundi]
    missed_quasicircular_distance = quasicircular_search_dist[quasicircular_search_missedi]

    emax = max(emax, max(quasicircular_search_e0))
    emin = min(emin, min(quasicircular_search_e0))

vt_ratios = []
vt_error_ratios = []
ifar_plotting = []
vt_values_qc = []
vt_values_ecc = []
vt_errors_qc = []
vt_errors_ecc = []

n_ecc_inj_in_coinc_time_ecc_search = len(eccentric_search_exclusive_ifar)
n_ecc_inj_in_coinc_time_qc_search = len(quasicircular_search_exclusive_ifar)

f_rec_inj_ecc_search = []
f_rec_inj_qc_search = []
f_rec_inj_ecc_search_err = []
f_rec_inj_qc_search_err = []

for ifar_threshold in ifar_thresholds:

    loud_ecc = eccentric_search_exclusive_ifar >= ifar_threshold
    quiet_ecc = eccentric_search_exclusive_ifar < ifar_threshold

    loud_qc = quasicircular_search_exclusive_ifar >= ifar_threshold
    quiet_qc = quasicircular_search_exclusive_ifar < ifar_threshold

    ecc_dist_found = found_eccentric_distance[loud_ecc]
    ecc_dist_missed = np.concatenate((found_eccentric_distance[quiet_ecc], missed_eccentric_distance))

    qc_dist_found = found_quasicircular_distance[loud_qc]
    qc_dist_missed = np.concatenate((found_quasicircular_distance[quiet_qc], missed_quasicircular_distance))

    v_total_ecc, v_total_ecc_err = sensitivity.volume_binned_pylal( ecc_dist_found, ecc_dist_missed,
                                                                        bins=vt_integration_bin )
    v_total_qc, v_total_qc_err = sensitivity.volume_binned_pylal( qc_dist_found, qc_dist_missed,
                                                                    bins=vt_integration_bin )

    vt_ratio_value = v_total_ecc/v_total_qc
    # foreground time is the same for both searches, 
    # so we can directly take the ratio of the sensitive volumes to get the ratio of the sensitive spacetime volumes. 
    # We only plot the ratio and its error if the ratio is finite.

    if vt_ratio_value !=float('inf'):
        vr_error = vratio_error(v_total_ecc, v_total_qc, v_total_ecc_err, v_total_qc_err)
        vt_ratios.append(v_total_ecc/v_total_qc)
        vt_error_ratios.append(vr_error)
        vt_values_qc.append(v_total_qc)
        vt_values_ecc.append(v_total_ecc)
        vt_errors_qc.append(v_total_qc_err)
        vt_errors_ecc.append(v_total_ecc_err)
        ifar_plotting.append(ifar_threshold)
        f_rec_inj_ecc_search.append(sum(loud_ecc)/n_ecc_inj_in_coinc_time_ecc_search)
        f_rec_inj_qc_search.append(sum(loud_qc)/n_ecc_inj_in_coinc_time_qc_search)
        f_rec_inj_ecc_search_err.append(np.sqrt(sum(loud_ecc))/n_ecc_inj_in_coinc_time_ecc_search)
        f_rec_inj_qc_search_err.append(np.sqrt(sum(loud_qc))/n_ecc_inj_in_coinc_time_qc_search)


fig, ax1  = plt.subplots()
eb = ax1.errorbar(ifar_plotting, vt_ratios, yerr=vt_error_ratios, color='magenta')#, fmt='o')
ax1.set_xlabel("IFAR [yr]")
ax1.set_ylabel(r"$\frac{\langle VT \rangle_{\rm EC}^{}}{\langle VT \rangle_{\rm QC}^{}}$")
ax1.set_xscale('log')
ax2 = ax1.twinx()
l1, = ax2.plot(ifar_plotting, f_rec_inj_ecc_search,
        #yerr=f_rec_inj_ecc_search_err ,
        color='blue', ls= '-.')
l2,  = ax2.plot(ifar_plotting, f_rec_inj_qc_search,
        #yerr=f_rec_inj_qc_search_err,
        color='orange', ls= '--')
ax2.set_ylabel(r"$f_{\rm{recovered}}^{{\rm EC},{\rm QC}}$")
ax2.set_yscale('log')
ax2.yaxis.set_minor_formatter(matplotlib.ticker.FuncFormatter(log_minor_fmt))
handles = [eb.lines[0], l1, l2]
labels = [r'$\langle VT \rangle^{}_{\rm EC} / \langle VT \rangle^{}_{\rm QC}$',
    r'$f_{\rm recovered}^{\rm EC}$',
    r'$f_{\rm recovered}^{\rm QC}$']

fig.legend(
    handles,
    labels,
    loc='center left',
    bbox_to_anchor=(.15, .55),
    frameon=False,
)

fig.savefig("fig4.png", bbox_inches='tight')
fig.savefig("fig4.pdf", bbox_inches='tight')