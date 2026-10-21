import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from pycbc.io import HFile, SingleDetTriggers
import pycbc.conversions

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

# Hardcoded parameters and file locations
single_trigger_file = 'L1_triggers_eccentric_search.hdf'
ifo = 'L1'
bank_file = 'Eccentric_template_bank_GW200105_search.hdf'
statmap_file = 'H1L1V1-COMBINE_STATMAP_ECCENTRIC_SEARCH.hdf'
time_window = 1
gps_start_time = 1262276674
gps_end_time = 1262276694

# Extract clustered triggers around GW200105 from statmap file
f_statmap = HFile(statmap_file, 'r')
d_statmap = f_statmap['foreground']
stat_arg_max = np.argmax(d_statmap['L1/stat'][:])
gw200105_trigger_id = d_statmap['L1/L1/trigger_id'][:][stat_arg_max]
gw200105_trigger_time = d_statmap['L1/L1/time'][:][stat_arg_max]
f_statmap.close()

# Filter triggers within time window
with HFile(single_trigger_file, 'r') as data_file:
    idx, _ = data_file.select( lambda endtime: abs(endtime - gw200105_trigger_time) < time_window, 
                          'end_time', group = ifo, return_data=False  )
    data_mask = np.zeros(data_file[ifo]['snr'].size, dtype=bool)
    data_mask[idx] =  True

single_triggers = SingleDetTriggers( single_trigger_file, ifo,
                                        bank_file = bank_file, premask=data_mask)

# Store data in dictionary instead of HDF5 file
trigger_data = {
    'mass1': single_triggers.mass1,
    'mass2': single_triggers.mass2,
    'eccentricity': single_triggers.eccentricity,
    'snr': single_triggers.snr,
    'chisq': single_triggers.rchisq,
    'newsnr_sgveto_psdvar_threshold': single_triggers.newsnr_sgveto_psdvar_threshold
}

colormap = 'plasma'

eccentricity = trigger_data['eccentricity']
snr = trigger_data['snr']
chisq = trigger_data['chisq']
mass1 = trigger_data['mass1']
mass2 = trigger_data['mass2']
newsnr_sgveto_psdvar_threshold = trigger_data['newsnr_sgveto_psdvar_threshold']

eccentric_chirp_mass = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(mass1, mass2, eccentricity)

# Best matched template parameter from Table II of Phukon et al., arxiv:2512.10803
search_m1 = 9.18
search_m2 = 2.06
search_s1z = -0.08
search_s2z = 0.05
search_ecc = 0.105


search_mch_ecc = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(search_m1, search_m2, search_ecc)
search_chi_eff = pycbc.conversions.chi_eff(search_m1, search_m2, search_s1z, search_s2z)

snr_mask = snr>9
snr_chisq_mask = np.logical_and(snr>8, chisq<1.5)


idx = np.argsort(newsnr_sgveto_psdvar_threshold[snr_chisq_mask])

plt.figure()

scatter = plt.scatter(eccentricity[snr_chisq_mask][idx], eccentric_chirp_mass[snr_chisq_mask][idx], 
                    c = newsnr_sgveto_psdvar_threshold[snr_chisq_mask][idx], s=5, cmap=colormap)
cb = plt.colorbar(scatter, pad=0.05)
cb.set_label(r"$\hat{\rho}$")#, rotation=370)#, labelpad=0.05)
cb.mappable.set_clim(vmin=8, vmax=13 )

plt.scatter( search_ecc, search_mch_ecc, marker='x', color='black')

plt.xlabel(r"$e_{20}^{}$")
plt.ylabel(r"$\mathcal{M}_c^{\rm ecc}$")

plt.xscale('log')

plt.tight_layout()
plt.savefig("fig7.png")
plt.savefig("fig7.pdf")