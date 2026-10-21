import numpy, matplotlib
from matplotlib import colors
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from pycbc.io import (
    get_chisq_from_file_choice, SingleDetTriggers, HFile
)
from pycbc import conversions
from pycbc.detector import Detector
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.collections import PathCollection
import math

rc_params = {
    # Color scheme (from Mathematica)
    'axes.prop_cycle': plt.cycler('color', ['5d81b4', 'e09b24', '8eb031', 'eb6235', '8678b2', 'c46e1a', '5c9dc7', 'ffbf00', 'a5609c']) + plt.cycler('ls', ['-', '--', '-.', (0, (1,0.85)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), (0, (5, 1)), ':', (4, (10, 3))]),
    
    # Axes
    'axes.titlesize': 24,
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
    'xtick.labelsize': 18,
    
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
    'ytick.labelsize': 18,
    
    # Legend
    'legend.frameon': True,
    'legend.fontsize': 12,
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
    'font.serif': ['Times New Roman'],
    'font.sans-serif': ['Bitstream Vera Sans'],
    'font.size': 24,
    'text.latex.preamble': r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts} \usepackage{xcolor}',
    
    # Saving figures
    'path.simplify': True,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'savefig.dpi': 200,
    
    # Additional customizations
    'backend': 'pdf',
    'patch.force_edgecolor': True,
}

plt.rcParams.update(rc_params)

def snr_from_chisq(chisq, newsnr, q=6.):
    snr = numpy.zeros(len(chisq)) + float(newsnr)
    ind = numpy.where(chisq > 1.)[0]
    snr[ind] = float(newsnr) / ( 0.5 * (1. + chisq[ind] ** (q/2.)) ) ** (-1./q)
    return snr

def bkg_snr_chisq(trigs, b_tids, chisq_choice='traditional'):
    bkg_snr = trigs.snr
    bkg_chisq = get_chisq_from_file_choice(trigs, chisq_choice)

    bkg_pos = numpy.logical_and(bkg_chisq > 0,  bkg_snr >=6)
    bkg_snr = bkg_snr[bkg_pos]
    bkg_chisq = bkg_chisq[bkg_pos]
    
    return bkg_snr, bkg_chisq

def trigs_snr_chisq(trigs, tids, chisq_choice='traditional'):
    snr = trigs.snr
    chisq = get_chisq_from_file_choice(trigs, chisq_choice=chisq_choice)

    _, tid_idx_inverse = numpy.unique(tids, return_inverse=True)

    snr = snr[tid_idx_inverse]
    chisq = chisq[tid_idx_inverse]
    pos = numpy.logical_and(chisq > 0, snr >= 6)
    
    return snr[pos], chisq[pos], pos

# Function to build a black proxy handle for a plotted artist
def make_black_proxy(handle, marker=None, markersize=None):
    if isinstance(handle, Line2D):
        # preserve marker shape and linestyle but force black color
        marker = handle.get_marker()
        ls = handle.get_linestyle()
        ms = handle.get_markersize() or 6
        lw = handle.get_linewidth() or 1.0
        # For pure marker-only lines, linestyle may be 'None' -> keep that
        return Line2D([0], [0],
                      color='black',
                      marker=marker if marker not in (None, 'None') else 'o',
                      linestyle=ls if ls not in (None, 'None') else 'None',
                      markersize=ms,
                      linewidth=lw)
    if isinstance(handle, PathCollection):
        # scatter: use a marker proxy. 's' in scatter are areas (points^2), so sqrt to get marker size
        sizes = handle.get_sizes()
        size = sizes[0] if len(sizes) else 36  # default area
        if markersize is None:
            markersize = math.sqrt(size)
        # try to infer edge width/marker shape if needed; here we use a filled circle
        return Line2D([0], [0],
                      color='black',
                      marker=marker if marker not in (None, 'None') else 'o',
                      linestyle='None',
                      markersize=markersize,
                      markeredgecolor='black')
    # fallback proxy
    return Line2D([0], [0], color='black', marker='o', linestyle='None')

# Statistics files of clustered background triggers in the eccentric search
h1v1_ecc_coinc_statistic_file = 'H1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf'
h1l1_ecc_coinc_statistic_file = 'H1L1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf'
l1v1_ecc_coinc_statistic_file = 'L1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_ECCENTRIC_SEARCH.hdf' 
h1l1v1_ecc_coinc_statistic_file = 'H1L1V1-EXCLUDE_ZEROLAG_FULL_DATA_3DET_ECCENTRIC_SEARCH.hdf' 

# Statistics files of clustered background triggers in the quasicircular search
h1v1_qc_coinc_statistic_file = 'H1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf'
h1l1_qc_coinc_statistic_file = 'H1L1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf'
l1v1_qc_coinc_statistic_file = 'L1V1-EXCLUDE_ZEROLAG_FULL_DATA_2DET_QUASICIRCULAR_SEARCH.hdf' 
h1l1v1_qc_coinc_statistic_file = 'H1L1V1-EXCLUDE_ZEROLAG_FULL_DATA_3DET_QUASICIRCULAR_SEARCH.hdf' 

ecc_stat_files_dict = {
                    "H1":[h1v1_ecc_coinc_statistic_file, h1l1_ecc_coinc_statistic_file, h1l1v1_ecc_coinc_statistic_file],
                    "L1":[h1l1_ecc_coinc_statistic_file, l1v1_ecc_coinc_statistic_file, h1l1v1_ecc_coinc_statistic_file],
                    "V1":[h1v1_ecc_coinc_statistic_file, l1v1_ecc_coinc_statistic_file, h1l1v1_ecc_coinc_statistic_file]
                    }

qc_stat_files_dict = {
                    "H1":[h1v1_qc_coinc_statistic_file, h1l1_qc_coinc_statistic_file, h1l1v1_qc_coinc_statistic_file],
                    "L1":[h1l1_qc_coinc_statistic_file, l1v1_qc_coinc_statistic_file, h1l1v1_qc_coinc_statistic_file],
                    "V1":[h1v1_qc_coinc_statistic_file, l1v1_qc_coinc_statistic_file, h1l1v1_qc_coinc_statistic_file]
                    }

# Single detector trigger files from the eccentric search
# Files are too large to be included in the data release
h1_single_ecc_trigger_file = 'H1_triggers_with_background_eccentric_search.hdf'
l1_single_ecc_trigger_file = 'L1_triggers_with_background_eccentric_search.hdf'
v1_single_ecc_trigger_file = 'V1_triggers_with_background_eccentric_search.hdf'

eccentric_single_trigger_files = {
                    "H1": h1_single_ecc_trigger_file,
                    "L1": l1_single_ecc_trigger_file,
                    "V1": v1_single_ecc_trigger_file
                    }

eccentric_bkg_tigger_ids = {
        "H1": numpy.array([]),
        "L1": numpy.array([]),
        "V1": numpy.array([])
        }

# Single detector trigger files with background from the quasicircular search
#  Files are too large to be included in the data release 
h1_single_qc_trigger_file = 'H1_triggers_with_background_quasicircular_search.hdf' 
l1_single_qc_trigger_file = 'L1_triggers_with_background_quasicircular_search.hdf'  
v1_single_qc_trigger_file = 'V1_triggers_with_background_quasicircular_search.hdf' 

qc_single_trigger_files = {
                    "H1": h1_single_qc_trigger_file,
                    "L1": l1_single_qc_trigger_file,
                    "V1": v1_single_qc_trigger_file
                    }

qc_bkg_tigger_ids = {
        "H1": numpy.array([]),
        "L1": numpy.array([]),
        "V1": numpy.array([])
        }

# Trigger files for quasi-circular injections in the eccentric search and quasicircular search
eccentric_single_injection_files = {
                    "H1": H1-QUASICIRCULAR_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hdf,
                    "L1": L1-QUASICIRCULAR_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hdf,
                    "V1": V1-QUASICIRCULAR_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hdf
        }

qc_single_injection_files = {
                    "H1": H1-QUASICIRCULAR_INJECTION_TRIGGERS_QUASICIRCULAR_SEARCH.hdf,
                    "L1": L1-QUASICIRCULAR_INJECTION_TRIGGERS_QUASICIRCULAR_SEARCH.hdf,
                    "V1": V1-QUASICIRCULAR_INJECTION_TRIGGERS_QUASICIRCULAR_SEARCH.hdf
        }

alpha = .9

for i, ifo in enumerate(['H1', 'L1', 'V1']):
    coinc_list = ecc_stat_files_dict[ifo]
    for j, coinc_statistic_file in enumerate(coinc_list):
        with HFile(coinc_statistic_file, 'r') as csf:
            b_tids = csf['background_exc'][ifo]['trigger_id'][:]
            b_tids = b_tids[b_tids >= 0]
            ifos = csf.attrs['ifos'].split(' ')
            eccentric_bkg_tigger_ids[ifo] = numpy.append(eccentric_bkg_tigger_ids[ifo].astype('i'), b_tids)
    
    coinc_list = qc_stat_files_dict[ifo]
    for j, coinc_statistic_file in enumerate(coinc_list):
        with HFile(coinc_statistic_file, 'r') as csf:
            b_tids = csf['background_exc'][ifo]['trigger_id'][:]
            b_tids = b_tids[b_tids >= 0]
            ifos = csf.attrs['ifos'].split(' ')
            qc_bkg_tigger_ids[ifo] = numpy.append(qc_bkg_tigger_ids[ifo].astype('i'), b_tids)

for key in eccentric_bkg_tigger_ids:
    eccentric_bkg_tigger_ids[key] = numpy.unique(eccentric_bkg_tigger_ids[key])

for key in qc_bkg_tigger_ids:
    qc_bkg_tigger_ids[key] = numpy.unique(qc_bkg_tigger_ids[key])

total_qc_trigger = 0
total_ec_trigger = 0

for ifo in ["H1", "L1", "V1"]:
    print (f"eccentric_bkg_tigger_ids {ifo}, {len(eccentric_bkg_tigger_ids[ifo])}")
    total_ec_trigger += len(eccentric_bkg_tigger_ids[ifo])
    print (f"qc_bkg_tigger_ids {ifo}, {len(qc_bkg_tigger_ids[ifo])}")
    total_qc_trigger += len(qc_bkg_tigger_ids[ifo])

scatter_size = 10
scatter_size_inj = 10

colormap = 'plasma'

cm = 1/2.54
fig, ax = plt.subplots(1,2, figsize=(23.2*cm, 13.5*cm), sharey=True)
fig.subplots_adjust(wspace=0.02)

# Add the eccentric background triggers

snr_min = []
chisq_min = []

snr_max = []
chisq_max = []

ecc_inj_snrs = numpy.array([])
ecc_inj_chisq = numpy.array([])
ecc_inj_e0s = numpy.array([])
ecc_inj_poss = numpy.array([])

qc_inj_snrs = numpy.array([])
qc_inj_chisq = numpy.array([])
qc_inj_e0s = numpy.array([])
qc_inj_poss = numpy.array([])


ecc_snrs = numpy.array([])
ecc_chisqs = numpy.array([])

qc_snrs = numpy.array([])
qc_chisqs = numpy.array([])

# Summary of quasicircular injections in the eccentric search
found_ecc_injection_file = 'H1L1V1_INJ_FIND_QUASICIRCULAR_INJECTIONS_ECCENTRIC_SEARCH.hdf'

# Summary of quasicircular injections in the quasicircular search
found_qc_injection_file = 'H1L1V1_INJ_FIND_QUASICIRCULAR_INJECTIONS_QUASICIRCULAR_SEARCH.hdf'

for i, ifo in enumerate(['H1', 'L1', 'V1']):

    ecc_trigs_bkg = SingleDetTriggers( eccentric_single_trigger_files[ifo],
                                        ifo, premask=eccentric_bkg_tigger_ids[ifo]
                                    )

    ecc_snr_bkg, ecc_chisq_bkg =  bkg_snr_chisq(ecc_trigs_bkg, eccentric_bkg_tigger_ids[ifo])
    ecc_snrs = numpy.append(ecc_snrs, ecc_snr_bkg)
    ecc_chisqs = numpy.append(ecc_chisqs, ecc_chisq_bkg)

    del ecc_trigs_bkg
    del ecc_snr_bkg
    del ecc_chisq_bkg


    qc_trigs_bkg = SingleDetTriggers( qc_single_trigger_files[ifo],
                                        ifo, premask=qc_bkg_tigger_ids[ifo]
                                    )

    qc_snr_bkg, qc_chisq_bkg =  bkg_snr_chisq(qc_trigs_bkg, qc_bkg_tigger_ids[ifo])
    qc_snrs = numpy.append(qc_snrs, qc_snr_bkg)
    qc_chisqs = numpy.append(qc_chisqs, qc_chisq_bkg)

    del qc_trigs_bkg
    del qc_snr_bkg
    del qc_chisq_bkg


    f = HFile(found_ecc_injection_file, 'r')
    inj_tid = f['found_after_vetoes'][ifo]['trigger_id'][:]
    valid_inj = inj_tid >= 0
    inj_tid = inj_tid[valid_inj]

    eff_dists = Detector(ifo).effective_distance(f['injections/distance'][:],
                                             f['injections/ra'][:],
                                             f['injections/dec'][:],
                                             f['injections/polarization'][:],
                                             f['injections/tc'][:],
                                             f['injections/inclination'][:])

    inj_idx = f['found_after_vetoes/injection_index'][valid_inj]
    eff_dist = eff_dists[inj_idx]
    m1, m2 = f['injections/mass1'][:][inj_idx], f['injections/mass2'][:][inj_idx]
    s1, s2 = f['injections/spin1z'][:][inj_idx], f['injections/spin2z'][:][inj_idx]
    
    mchirp = conversions.mchirp_from_mass1_mass2(m1, m2)
    eta = conversions.eta_from_mass1_mass2(m1, m2)
    weighted_spin = conversions.chi_eff(m1, m2, s1, s2)

    if len(inj_tid):
        inj_trigs = SingleDetTriggers(
            eccentric_single_injection_files[ifo],
            ifo, premask=inj_tid
        )
        inj_snr, inj_chisq, inj_pos = trigs_snr_chisq(inj_trigs, inj_tid)

        del inj_trigs
    else:
        inj_snr = numpy.array([])
        inj_chisq = numpy.array([])
        inj_pos = numpy.array([])

    ecc_inj_snrs = numpy.append(ecc_inj_snrs, inj_snr)
    ecc_inj_chisq = numpy.append(ecc_inj_chisq, inj_chisq)
    ecc_inj_poss = numpy.append(ecc_inj_poss, inj_pos)
 
    f = HFile(found_qc_injection_file, 'r')
    inj_tid = f['found_after_vetoes'][ifo]['trigger_id'][:]
    valid_inj = inj_tid >= 0
    inj_tid = inj_tid[valid_inj]

    eff_dists = Detector(ifo).effective_distance(f['injections/distance'][:],
                                             f['injections/ra'][:],
                                             f['injections/dec'][:],
                                             f['injections/polarization'][:],
                                             f['injections/tc'][:],
                                             f['injections/inclination'][:])

    inj_idx = f['found_after_vetoes/injection_index'][valid_inj]
    eff_dist = eff_dists[inj_idx]
    m1, m2 = f['injections/mass1'][:][inj_idx], f['injections/mass2'][:][inj_idx]
    s1, s2 = f['injections/spin1z'][:][inj_idx], f['injections/spin2z'][:][inj_idx]

    mchirp = conversions.mchirp_from_mass1_mass2(m1, m2)
    eta = conversions.eta_from_mass1_mass2(m1, m2)
    weighted_spin = conversions.chi_eff(m1, m2, s1, s2)

    if len(inj_tid):
        inj_trigs = SingleDetTriggers(
            qc_single_injection_files[ifo],
            ifo, premask=inj_tid
        )
        inj_snr, inj_chisq, inj_pos = trigs_snr_chisq(inj_trigs, inj_tid)
        
        del inj_trigs

    else:
        inj_snr = numpy.array([])
        inj_chisq = numpy.array([])
        inj_pos = numpy.array([])

    qc_inj_snrs = numpy.append(qc_inj_snrs, inj_snr)
    qc_inj_chisq = numpy.append(qc_inj_chisq, inj_chisq)
    qc_inj_poss = numpy.append(qc_inj_poss, inj_pos)



bkg_qc_line = ax[0].scatter(qc_snrs, qc_chisqs, marker='x', color='black',
                            label='Background: QC', s = scatter_size_inj)
bkg_ec_line = ax[1].scatter(ecc_snrs, ecc_chisqs, marker='x', color='black',
                            label='Background: EC', s = scatter_size_inj)

inj_line1 = ax[1].scatter(ecc_inj_snrs, ecc_inj_chisq,
                  s = scatter_size_inj,#c=ecc_inj_e0s,
                  marker='^', linewidth=0, label="Injections: EC",
                  cmap=colormap, alpha=alpha)
inj_line2 = ax[0].scatter(qc_inj_snrs, qc_inj_chisq, #c=qc_inj_e0s,
                  s=scatter_size_inj,
                  marker='^', linewidth=0, label="Injections: QC",
                  cmap=colormap, alpha=alpha)

inj_line2_proxy = make_black_proxy(inj_line2, marker='^', markersize=numpy.sqrt(scatter_size_inj))

# SNR and chi-squared values for the GW200105 trigger in the eccentric search
gw200105_line1 = ax[1].scatter(13.41, 1.14, 
        marker='D', color='red', s=2*scatter_size_inj)

# SNR and chi-squared values for the GW200105 trigger in the quasi-circular search        
gw200105_line2 = ax[0].scatter(13.42, 1.2, 
        marker='D', color='red',  s=2*scatter_size_inj) 

gw200105_line2_proxy = make_black_proxy(gw200105_line2, marker='D', markersize=numpy.sqrt(scatter_size_inj))

ax[0].set_ylabel(r"$\chi_r^2$")
ax[0].set_xlabel(r"SNR, $\rho$")
ax[1].set_xlabel(r"SNR, $\rho$")

ax[1].set_xlim([6,100])

ax[0].set_xlim([6,100])
ax[0].set_ylim([1,100])

r =  numpy.logspace( 1, 100, 200)

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_xscale('log')
ax[1].set_xscale('log')

lines_qc = [(bkg_qc_line), (inj_line2_proxy), (gw200105_line2_proxy)]

labels_qc = ["Noise","Injections",r"GW200105"]

lines_ec = [(bkg_ec_line) ,(inj_line1),(gw200105_line1)]
labels_ec = ["Noise","Injections",r"GW200105"]

ax[0].set_title("Quasi-circular search")
ax[1].set_title("Eccentric search")

legend1 = ax[0].legend(lines_qc, labels_qc,
            ncol=1, loc='upper left',
            handler_map={tuple: HandlerTuple(ndivide=None)},
            framealpha=0.5, frameon=False,
            edgecolor='black',
            markerfirst=True,
            markerscale=2)


for handle in legend1.legendHandles:
    handle.set_color('black')

for text in legend1.get_texts():
    text.set_color('black')


fig.savefig('fig3.png', bbox_inches="tight")
fig.savefig('fig3.pdf', bbox_inches="tight")

