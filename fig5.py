import matplotlib as mpl; mpl.use('Agg')
import numpy as np
import h5py
from matplotlib import pyplot as plt
from scipy.stats import norm, poisson
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pycbc.results
from pycbc import conversions
from pycbc.io.hdf import HFile
from pycbc.events import significance
from matplotlib import rcParams

rc_params = {
    # Color scheme (from Mathematica)
    'axes.prop_cycle': plt.cycler('color', ['5d81b4', 'e09b24', '8eb031', 'eb6235', '8678b2', 'c46e1a', '5c9dc7', 'ffbf00', 'a5609c']) + plt.cycler('ls', ['-', '--', '-.', (0, (1,0.85)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), (0, (5, 1)), ':', (4, (10, 3))]),
    
    # Axes
    'axes.titlesize': 28,
    'axes.labelsize': 24,
    'axes.formatter.use_mathtext': True,
    'axes.linewidth': 1.75,
    
    # Grid lines
    'grid.linewidth': 0.5,
    'grid.linestyle': 'dashed',
    'grid.color': 'xkcd:light gray',
    
    # Lines
    'lines.linewidth': 1,
    'lines.markersize': 10,
    'hatch.linewidth': 0.25,
    'patch.antialiased': True,
    
    # Ticks
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
    'xtick.labelsize': 10,
    
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
    'ytick.labelsize': 10,
    
    # Legend
    'legend.frameon': True,
    'legend.fontsize': 18,
    'legend.handlelength': 1.375,
    'legend.labelspacing': 0.4,
    'legend.columnspacing': 1,
    'legend.facecolor': 'white',
    'legend.edgecolor': 'white',
    'legend.framealpha': 1,
    'legend.title_fontsize': 8,
    
    # Figure size
    'figure.figsize': (3.25, 2.0),
    'figure.subplot.left': 0.125,
    'figure.subplot.bottom': 0.175,
    'figure.subplot.top': 0.95,
    'figure.subplot.right': 0.95,
    'figure.autolayout': False,
    
    # Fonts
    'text.usetex': True,
    'mathtext.fontset': 'stixsans',
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman', 'cmr10'],
    'font.sans-serif': ['Bitstream Vera Sans'],
    'font.size': 10,
    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts} \usepackage{xcolor}',
    
    # Saving figures
    'path.simplify': True,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'savefig.dpi': 400,
    
    # Additional customizations
    'backend': 'pdf',
    'patch.force_edgecolor': True,
}

plt.rcParams.update(rc_params)

# Clustered trigger files
eccentric_trigger_files = [
    'H1L1V1-COMBINE_STATMAP_ECCENTRIC_SEARCH.hdf'
]
quasicircular_trigger_files = [
    'H1L1V1-COMBINE_STATMAP_QUASICIRCULAR_SEARCH.hdf'
]

eccentric_statmap_files = [
    'L1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf',
    'H1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf',
    'H1L1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf',
    'H1L1V1-EXCLUDE_ZEROLAG_FULL_DATA_3DET_ECCENTRIC_SEARCH.hdf',
    'V1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_ECCENTRIC_SEARCH.hdf',
    'L1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_ECCENTRIC_SEARCH.hdf',
    'H1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_ECCENTRIC_SEARCH.hdf',
]

quasicircular_statmap_files = [
    'L1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf',
    'H1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf',
    'H1L1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf',
    'H1L1V1-EXCLUDE_ZEROLAG_FULL_DATA_3DET_QUASICIRCULAR_SEARCH.hdf',
    'V1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_QUASICIRCULAR_SEARCH.hdf',
    'L1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_QUASICIRCULAR_SEARCH.hdf',
    'H1-EXCLUDE_ZEROLAG_FULL_DATA_1DET_QUASICIRCULAR_SEARCH.hdf',
]

min_x = -2.0
max_x = 30.0

ifo_combos = ['L1', 'L1V1', 'H1L1', 'H1V1', 'L1V1', 'H1L1V1']

# Significance options (for trigger_fit method)
significance_options = {
    'H1': {
        'limit_ifar': 1000,
        'far_calculation_method': 'trigger_fit',
        'fit_function': 'exponential',
        'fit_threshold': -5.5,
    },
    'L1': {
        'limit_ifar': 1000,
        'far_calculation_method': 'trigger_fit',
        'fit_function': 'exponential',
        'fit_threshold': -4.3,
    },
    'V1': {
        'limit_ifar': 1000,
        'far_calculation_method': 'trigger_fit',
        'fit_function': 'exponential',
        'fit_threshold': -8.3,
    },
}

coinc_color_map = {
    'H':  '#ee0000',  # red
    'L':  'blue',#'#4ba6ff',  # blue
    'V':  '#9b59b6',  # magenta/purple
    'HL': '#0000DD',
    'LV': 'black',#'#DFFF00',#'#7FFF00', #chartreuse
    'HV': '#990000',
    'HLV': np.array([20, 120, 20]) / 255.}

coinc_color_map_qc = {
    'H':  '#ee0000',  # red
    'L':  '#691CE3', #'#96E31C',#'orange',
    'V':  '#9b59b6',  # magenta/purple
    'HL': '#0000DD',
    'LV': 'magenta',
    'HV': '#990000',
    'HLV': np.array([20, 120, 20]) / 255.}

coinc_line_map = {
    'H':  '--',
    'L':  ':',
    'V':  '--',
    'HL': '-',
    'LV': 'dashed',
    'HV': '-',
    'HLV':'-'
}

coinc_line_map_qc = {
    'H':   ':',
    'L':  (0, (3, 5, 1, 5, 1, 5)),# (10, (3, 1, 1, 1, 1, 1)),#(10, (3, 10, 1, 10, 1, 10)),
    'V':   ':',
    'HL':  '-.',
    'LV':  '-.',
    'HV':  '-.',
    'HLV': '-.'
}


ecc_foreground_color = 'red'
qc_foreground_color = 'black'

# Marker size for foreground triggers
ms=5

trigf_ecc = [HFile(f, 'r') for f in eccentric_trigger_files]
trigf_qc = [HFile(f, 'r') for f in quasicircular_trigger_files]

# Using exclusive IFAR values for plotting and cumulative number calculation
ifar_str = 'ifar_exc'

# get foreground IFAR values and cumulative number for each IFAR value
fore_ifar_ecc = np.array([])
for f in trigf_ecc:
    fore_ifar_ecc = np.append(fore_ifar_ecc, f['foreground/' + ifar_str][:])
    
fore_ifar_qc = np.array([])
for f in trigf_qc:
    fore_ifar_qc = np.append(fore_ifar_qc, f['foreground/' + ifar_str][:])


sortifar_ecc = fore_ifar_ecc.argsort()
fore_ifar_ecc = fore_ifar_ecc[sortifar_ecc]
fore_cumnum_ecc = np.arange(len(fore_ifar_ecc), 0, -1)

sortifar_qc = fore_ifar_qc.argsort()
fore_ifar_qc = fore_ifar_qc[sortifar_qc]
fore_cumnum_qc = np.arange(len(fore_ifar_qc), 0, -1)

# get expected (noise-only) foreground IFAR values and cumulative number for
# each IFAR value
expected_ifar = np.logspace(-4., np.log10(10000),
                               num=1000, base=10.0)
fg_time_ecc = 0
for f in trigf_ecc:
    fg_time_ecc += f.attrs['foreground_time']

expected_cumnum_ecc = conversions.sec_to_year(fg_time_ecc) / expected_ifar

fg_time_qc = 0
for f in trigf_qc:
    fg_time_qc += f.attrs['foreground_time']


print ("fg_time_ecc %s, fg_time_qc %s"%(fg_time_ecc, fg_time_qc))

expected_cumnum_qc = conversions.sec_to_year(fg_time_qc) / expected_ifar


# make figure
cm = 1/2.54
fig, ax = plt.subplots(1,2, figsize=(23.2*cm, 13.5*cm))
fig.subplots_adjust(wspace=0.3, top=0.91)

divider = make_axes_locatable(ax[0])
axLin = divider.append_axes("right", size=1.0, pad=0.05, sharey=ax[0])

# plot the expected background
ax[0].plot(expected_ifar, expected_cumnum_ecc, linestyle='solid', linewidth=2,
        color='black', label='Expected background')

lower = {}
upper = {}
for sigma in [1, 2, 3, 4]:
    # one-tailed significance
    eps = norm.sf(sigma, 0, 1)
    upper[sigma] = []
    lower[sigma] = []
    for mu in expected_cumnum_ecc[::-1]:
        # initialize at 'small' value
        nlo = max(0, int(mu - 5 * mu**0.5))
        while True:
            # cdf is sum of probabilities *including* nlo
            cdf = poisson.cdf(nlo, mu)
            if cdf > eps:
                break
            # count up
            nlo += 1
        # if value is not in tail then back off by one
        # nlo - 1 is largest value such that it is still in eps tail
        # thus put boundary at nlo - 0.5
        lower[sigma].append(nlo - 0.5)
        # initialize at 'large' value
        nup = max(10, int(mu + 5 * mu**0.5))
        while True:
            # sf is sum of probabilities *above* nup, ie complement of cdf
            surviv = poisson.sf(nup, mu)
            if surviv > eps:
                break
            # count down
            nup -= 1
        # here sum of probabilities from nup + 1 upwards is above eps
        # thus smallest value in eps tail is nup + 2
        # thus put boundary at nup + 1.5
        upper[sigma].append(nup + 1.5)

plotifar_ecc = expected_ifar[::-1]
ax[0].fill_between(plotifar_ecc, lower[4], lower[3], facecolor='k', alpha=0.15, edgecolor=None, linewidth=0)#,
                   #label=r'$<4\sigma$')
ax[0].fill_between(plotifar_ecc, lower[3], lower[2], facecolor='k', alpha=0.3, edgecolor=None, linewidth=0)#,
                   #label=r'$<3\sigma$')edgecolor=None
ax[0].fill_between(plotifar_ecc, lower[2], lower[1], facecolor='k', alpha=0.45, edgecolor=None, linewidth=0)#,
                   #label=r'$<2\sigma$')
ax[0].fill_between(plotifar_ecc, lower[1], upper[1], facecolor='k', alpha=0.6, edgecolor=None, linewidth=0)
                   #label=r'$<1\sigma$')
ax[0].fill_between(plotifar_ecc, upper[1], upper[2], facecolor='k', alpha=0.45, edgecolor=None, linewidth=0)
ax[0].fill_between(plotifar_ecc, upper[2], upper[3], facecolor='k', alpha=0.3, edgecolor=None, linewidth=0)
ax[0].fill_between(plotifar_ecc, upper[3], upper[4], facecolor='k', alpha=0.15, edgecolor=None, linewidth=0)

# plot the foreground triggers

fore_ifar_ecc_mask = fore_ifar_ecc<fore_ifar_ecc[-1]
ax[0].plot(fore_ifar_ecc[fore_ifar_ecc_mask], fore_cumnum_ecc[fore_ifar_ecc_mask], linestyle='None', color=ecc_foreground_color,
        marker='+', ms=ms, label='Foreground: EC')
axLin.plot(fore_ifar_ecc, fore_cumnum_ecc, linestyle='None', color=ecc_foreground_color,
        marker='+', ms=ms)
max_ifar_ecc = max(fore_ifar_ecc)

fore_ifar_qc_mask = fore_ifar_qc<fore_ifar_qc[-1]
ax[0].plot(fore_ifar_qc[fore_ifar_qc_mask], fore_cumnum_qc[fore_ifar_qc_mask], linestyle='None', color=qc_foreground_color,
        marker='x', ms=ms, label='Foreground: QC')
axLin.plot(fore_ifar_qc, fore_cumnum_qc, linestyle='None', color=qc_foreground_color,
        marker='x', ms=ms)
max_ifar_qc = max(fore_ifar_qc)


# format plot
# scale the plot around the foreground
ax[0].set_ylim(0.7, 1e3)#0.05 * len(fore_cumnum_ecc))
print ("fore_ifar_ecc",min(fore_ifar_ecc))
print ("max_ifar", max_ifar_ecc)
if min(fore_ifar_ecc) < 1e-5:
    ax[0].set_xlim(1e-5, 3 * max_ifar_ecc)
else:
    ax[0].set_xlim(30 * min(fore_ifar_ecc), 3 * max_ifar_ecc)

ax[0].set_ylabel('Event count')#, size='large')
ifar_label = 'IFAR [yr]'#'Inverse False Alarm Rate (yr)'
ax[0].set_xlabel(ifar_label, loc='right')

ax[0].set_yscale('log')
ax[0].set_xscale('log')
axLin.set_yscale('log')
ax[0].set_xlim(xmax=948)
axLin.set_xlim(xmin=948, xmax=1010)

ax[0].spines['right'].set_visible(False)
ax[0].yaxis.set_ticks_position('left')
axLin.spines['left'].set_visible(False)
axLin.yaxis.set_ticks_position('right')
plt.setp(axLin.get_yticklabels(), visible=False)


# Create significance dictionary from detector options
sig_dict = {}
for combo in ifo_combos:
    # Extract individual detectors from combo (e.g., 'H1L1' -> ['H1', 'L1'])
    detectors = combo.replace('1', '').split()
    if len(detectors) == 0:
        detectors = [combo[0]]  # Single detector case
    
    # Use options from the first detector in the combo
    primary_det = detectors[0]
    if primary_det in significance_options:
        sig_dict[combo] = {
            'method': 'trigger_fit',
            'far_calculation_method': significance_options[primary_det]['far_calculation_method'],
            'fit_threshold': significance_options[primary_det]['fit_threshold'],
            'fit_function': significance_options[primary_det]['fit_function'],
            'limit_ifar': significance_options[primary_det]['limit_ifar'],
        }
    else:
        # Default values if detector not found
        sig_dict[combo] = {
            'method': 'trigger_fit',
            'far_calculation_method': 'trigger_fit',
            'fit_threshold': -5.0,
            'fit_function': 'exponential',
            'limit_ifar': 1000,
        }

ifar_calc_points = np.linspace(min_x, max_x, 200)

fitted = []  # List of combos for which a fit is used

for f in eccentric_statmap_files:
    with h5py.File(f, 'r') as exc_zlag_f:
        ifos_list = exc_zlag_f.attrs['ifos'].split()
        sig_key = exc_zlag_f.attrs['ifos'].replace(' ', '')  # eg 'H1L1'
        coinc_key = sig_key.replace('1', '')  # eg 'HL'
        stat_exc = exc_zlag_f['background_exc']['stat'][:]
        stat = exc_zlag_f['background']['stat'][:]
        dec_fac_exc = exc_zlag_f['background_exc']['decimation_factor'][:]
        dec_fac = exc_zlag_f['background']['decimation_factor'][:]
        bg_time = conversions.sec_to_year(exc_zlag_f.attrs['background_time_exc'])
        fg_time = conversions.sec_to_year(exc_zlag_f.attrs['foreground_time'])
        ifar_file = exc_zlag_f['background_exc']['ifar'][:]
        cstat_fore = exc_zlag_f['foreground']['stat'][:]
        cstat_fore.sort()
        cstat_rate = np.arange(len(cstat_fore), 0, -1)/fg_time

    # If fitting / extrapolation is done
    if sig_dict[sig_key]['method'] == 'trigger_fit':
        fitted.append(sig_key)
        fit_thresh = sig_dict[sig_key]['fit_threshold']
        print ("sig_dict",sig_dict[sig_key])
        far_tuple = significance.get_far(stat_exc,
                                         ifar_calc_points,
                                         dec_fac_exc,
                                         bg_time,
                                         **sig_dict[sig_key])
        far = far_tuple[1]  # Version-proof: get_far returns 2 or 3 values

        if coinc_key in ['L', 'LV']:
            print (f"coinc_key {coinc_key}, coinc_line_map[coinc_key] {coinc_line_map[coinc_key]}")
            ax[1].plot(ifar_calc_points, far, linestyle=coinc_line_map[coinc_key],
                c=coinc_color_map[coinc_key], label=coinc_key+ ": EC", zorder=-5, linewidth=2)
        
        if coinc_key=='L':
            cstat_fore_ecc_gw200105 = np.array([], dtype=float)
            cstat_rate_ecc_gw200105 = np.array([], dtype=float)
            rm_idx = cstat_fore.argmax()
            cstat_fore_ecc_gw200105 = np.append(cstat_fore_ecc_gw200105, cstat_fore[rm_idx])
            cstat_rate_ecc_gw200105 = np.append(cstat_rate_ecc_gw200105, cstat_rate[rm_idx])

        del fit_thresh  # avoid variable hanging around

    else:
        far_tuple = significance.get_far(stat_exc,
                                      ifar_calc_points,
                                      dec_fac_exc,
                                      bg_time,
                                      method='n_louder')

        if coinc_key in ['L', 'LV']:
            ax[1].plot(ifar_calc_points, far_tuple[1], linestyle=coinc_line_map[coinc_key],
                c=coinc_color_map[coinc_key], label=coinc_key+ ": EC", zorder=-5, linewidth=2)

    del sig_key, coinc_key, far_tuple  # avoid variables hanging around


for f in quasicircular_statmap_files:
    with h5py.File(f, 'r') as exc_zlag_f:
        sig_key = exc_zlag_f.attrs['ifos'].replace(' ', '')
        coinc_key = sig_key.replace('1', '')
        stat_exc = exc_zlag_f['background_exc']['stat'][:]
        stat = exc_zlag_f['background']['stat'][:]
        print (f"Compare stat for {coinc_key} excl {max(stat_exc)}, full {max(stat)}" )
        dec_fac_exc = exc_zlag_f['background_exc']['decimation_factor'][:]
        dec_fac = exc_zlag_f['background']['decimation_factor'][:]
        bg_time = conversions.sec_to_year(exc_zlag_f.attrs['background_time_exc'])
        fg_time = conversions.sec_to_year(exc_zlag_f.attrs['foreground_time'])
        ifar_file = exc_zlag_f['background_exc']['ifar'][:]
        cstat_fore = exc_zlag_f['foreground']['stat'][:]
        cstat_fore.sort()
        cstat_rate = np.arange(len(cstat_fore), 0, -1)/fg_time
    
    # If fitting / extrapolation is done
    if sig_dict[sig_key]['method'] == 'trigger_fit':
        fit_thresh = sig_dict[sig_key]['fit_threshold']
        far_tuple = significance.get_far(stat_exc,
                                         ifar_calc_points,
                                         dec_fac_exc,
                                         bg_time,
                                         **sig_dict[sig_key])
        far = far_tuple[1]

        if coinc_key=='L':
            cstat_fore_gc_gw200105 = np.array([], dtype=float)
            cstat_rate_gc_gw200105 = np.array([], dtype=float)
            rm_idx = cstat_fore.argmax()
            cstat_fore_gc_gw200105 = np.append(cstat_fore_gc_gw200105, cstat_fore[rm_idx])
            cstat_rate_gc_gw200105 = np.append(cstat_rate_gc_gw200105, cstat_rate[rm_idx])

        if coinc_key in ['L', 'LV']:
            ax[1].plot(ifar_calc_points, far,
                linestyle=coinc_line_map_qc[coinc_key],
                c=coinc_color_map_qc[coinc_key], label=coinc_key + ": QC", zorder=-5, linewidth=2)
        del fit_thresh  # avoid variable hanging around

    else:
        far_tuple = significance.get_far(stat_exc,
                                      ifar_calc_points,
                                      dec_fac_exc,
                                      bg_time,
                                      method='n_louder')

        if coinc_key in ['L', 'LV']:
            ax[1].plot(ifar_calc_points, far_tuple[1], linestyle=coinc_line_map_qc[coinc_key],
                c=coinc_color_map_qc[coinc_key], label=coinc_key + ": QC", zorder=-5, linewidth=2)

    del sig_key, coinc_key, far_tuple  # avoid variables hanging around

ax[1].plot( cstat_fore_ecc_gw200105, cstat_rate_ecc_gw200105, linestyle='None', color=ecc_foreground_color,
            marker='+', ms=ms)
ax[1].plot( cstat_fore_gc_gw200105, cstat_rate_gc_gw200105, linestyle='None', color=qc_foreground_color,
            marker='x', ms=ms)

ax[1].semilogy()
ax[1].set_xlim([min_x, max_x])
ax[1].set_ylim([0.7, .5*1e6])
ax[1].tick_params(axis='both')
ax[1].set_xlabel(r'$\Lambda_s^{}$')
ax[1].set_ylabel('Event rate [yr$^{-1}$]')

lines_labels = [axs.get_legend_handles_labels() for axs in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
fig.legend(lines, labels,
            loc='upper center', ncol=3, handlelength=3, bbox_to_anchor=(0., 1.06, 1., .102))

fig.savefig("fig5.png", bbox_inches="tight")
fig.savefig("fig5.pdf", bbox_inches="tight")