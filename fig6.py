import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from pycbc.io.hdf import HFile
import pycbc.conversions
import corner
from pycbc.io import SingleDetTriggers

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

# Injection find summary file for eccentric injections in the eccentric search
injection_file = 'H1L1V1_INJ_FIND_ECCENTRIC_INJECTIONS_ECCENTRIC_SEARCH.hdf'
# Eccentric template bank file
bank_file = 'Eccentric_template_bank_GW200105_search.hdf'
single_trigger_files = {'H1': 'H1-ECCENTRIC_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hdf',
                        'L1': 'L1-ECCENTRIC_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hdf',
                        'V1': 'V1-ECCENTRIC_INJECTION_TRIGGERS_ECCENTRIC_SEARCH.hd'}

f = HFile(injection_file,'r')

inj = f['injections']

found_cols, found_names, found_formats = [], [], []

ifos = f.attrs['ifos'].split(' ')
print ("ifos", ifos)

found = f['found_after_vetoes']
idx = found['injection_index'][:]

outfile_data = {'inj': {}, 'stat': {}}

outfile_data['inj']['mass1'] = inj['mass1'][:][idx]
outfile_data['inj']['mass2'] = inj['mass2'][:][idx]
outfile_data['inj']['eccentricity'] = inj['eccentricity'][:][idx]
outfile_data['inj']['spin1z'] = inj['spin1z'][:][idx]
outfile_data['inj']['spin2z'] = inj['spin2z'][:][idx] 
outfile_data['inj']['optimal_snr_H1'] = inj['optimal_snr_H1'][:][idx]
outfile_data['inj']['optimal_snr_L1'] = inj['optimal_snr_L1'][:][idx]
outfile_data['inj']['optimal_snr_V1'] = inj['optimal_snr_V1'][:][idx]

outfile_data['stat']['ifar'] = found['ifar_exc'][:]
outfile_data['stat']['ranking_stat'] = found['stat'][:]

detectors = f.attrs['ifos'].split(' ')
keys = found.keys()

print ("Keys in found dataset", list(keys))

detectors_used = []

for det in detectors:
    if(det in keys):
        detectors_used.append(det)
        print (f"More checks in found[ifo].keys() {det} {found[det].keys()}")
        f_det_keys = list(found[det].keys())
        for key_v in f_det_keys:
            print (f"{det}, {key_v} {len(found[det][key_v][:])}")
            print (f"{key_v} {found[det][key_v][:][0]}")
print ("detectors_used",detectors_used)

ids = {detector:found[detector+'/trigger_id'][:] for detector in detectors_used}

for detector in ids:
    print (f"Checks -1 id in ids for {detector}")
    print ( f"numpy argwhere {np.argwhere(ids[detector]==-11)}" )
    print ( f"Testing first elements  {ids[detector]}")

if single_trigger_files:
    for ifo in single_trigger_files:
        outfile_data[ifo] = {}
        inj_triggers = SingleDetTriggers( single_trigger_files[ifo],
                                        ifo, #premask=ids[ifo],
                                        bank_file = bank_file
                                    )
        ids_ifo = np.array(ids[ifo])
        ids_na = np.argwhere(ids_ifo == -1)
        
        snr_vals = inj_triggers.snr[ids_ifo]
        snr_vals[ids_ifo == -1] = np.nan
        outfile_data[ifo]['snr'] = snr_vals
        
        chisq_vals = inj_triggers.rchisq[ids_ifo]
        chisq_vals[ids_ifo == -1] = np.nan
        outfile_data[ifo]['chisq'] = chisq_vals
        
        newsnr_vals = inj_triggers.newsnr[ids_ifo]
        newsnr_vals[ids_ifo == -1] = np.nan
        outfile_data[ifo]['newsnr'] = newsnr_vals

        newsnr_sgveto_psdvar_threshold_vals = inj_triggers.newsnr_sgveto_psdvar_threshold[ids_ifo]
        newsnr_sgveto_psdvar_threshold_vals[ids_ifo == -1] = np.nan
        outfile_data[ifo]['newsnr_sgveto_psdvar_threshold'] = newsnr_sgveto_psdvar_threshold_vals

        template_mass1_values = inj_triggers.mass1[ids_ifo]
        template_mass1_values[ids_ifo == -1] = np.nan

        if ifo=='H1':
            h1_m1 = template_mass1_values
        if ifo=='L1':
            l1_m1 = template_mass1_values
        if ifo=='V1':
            v1_m1 = template_mass1_values
        
        outfile_data[ifo]['mass1'] = template_mass1_values
        
        template_mass2_values = inj_triggers.mass2[ids_ifo]
        template_mass2_values[ids_ifo == -1] = np.nan
        outfile_data[ifo]['mass2'] = template_mass2_values

        template_ecc_values = inj_triggers.eccentricity[ids_ifo]
        template_ecc_values[ids_ifo == -1] = np.nan
        outfile_data[ifo]['eccentricity'] = template_ecc_values

        template_spin1z_values = inj_triggers.spin1z[ids_ifo]
        template_spin1z_values[ids_ifo == -1] = np.nan
        outfile_data[ifo]['spin1z'] = template_spin1z_values

        template_spin2z_values = inj_triggers.spin2z[ids_ifo]
        template_spin2z_values[ids_ifo == -1] = np.nan
        outfile_data[ifo]['spin2z'] = template_spin2z_values

f.close()

print ("np.allclose(h1_m1[~(np.isnan(h1_m1) | np.isnan(l1_m1))], l1_m1[~(np.isnan(h1_m1) | np.isnan(l1_m1))])",np.allclose(h1_m1[~(np.isnan(h1_m1) | np.isnan(l1_m1))], l1_m1[~(np.isnan(h1_m1) | np.isnan(l1_m1))]))
print ('np.allclose(h1_m1[~(np.isnan(h1_m1) | np.isnan(v1_m1))], v1_m1[~(np.isnan(h1_m1) | np.isnan(v1_m1))])',np.allclose(h1_m1[~(np.isnan(h1_m1) | np.isnan(v1_m1))], v1_m1[~(np.isnan(h1_m1) | np.isnan(v1_m1))]))
print ('np.allclose(l1_m1[~(np.isnan(l1_m1) | np.isnan(v1_m1))], v1_m1[~(np.isnan(l1_m1) | np.isnan(v1_m1))])',np.allclose(l1_m1[~(np.isnan(l1_m1) | np.isnan(v1_m1))], v1_m1[~(np.isnan(l1_m1) | np.isnan(v1_m1))]))

get_data = outfile_data

print (get_data.keys())
print (get_data['H1'].keys())
print (get_data['inj'].keys())
search_m1 = 9.18
search_m2 = 2.06
search_s1z = -0.08
search_s2z = 0.05
search_ecc = 0.105
search_mchirp = pycbc.conversions.mchirp_from_mass1_mass2(search_m1, search_m2)
search_chieff = pycbc.conversions.chi_eff(search_m1, search_m2, search_s1z, search_s2z)
search_mchirp_ecc = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(search_m1, search_m2, search_ecc)

# Paremeters of the median posterior value for GW200105 from G. Morras et al., Astrophys. J. Lett. 1000, L2 (2026)
mps_m1_source = 11.5
mps_m2_source = 1.5
mps_mchirp_source = 3.33
mps_ecc = 0.145
mps_z = 0.06
mps_chi_eff = 0.004
mps_mchirp_ecc_source = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(mps_m1_source, mps_m2_source, mps_ecc)
z_factor = 1 + mps_z

mps_m1 = mps_m1_source * z_factor
mps_m2 = mps_m2_source * z_factor
mps_mchirp = mps_mchirp_source * z_factor
mps_mchirp_ecc = mps_mchirp_ecc_source * z_factor


print (get_data.keys())
print (get_data['H1'].keys())
print (get_data['inj'].keys())

mass1_vstack = np.vstack([get_data['H1']['mass1'][:], get_data['L1']['mass1'][:], get_data['V1']['mass1'][:] ])
template_mass1 =  np.nanmax(mass1_vstack, axis=0)

mass2_vstack = np.vstack([get_data['H1']['mass2'][:], get_data['L1']['mass2'][:], get_data['V1']['mass2'][:] ])
template_mass2 =  np.nanmax(mass2_vstack, axis=0)

eccentricity_vstack = np.vstack([get_data['H1']['eccentricity'][:], get_data['L1']['eccentricity'][:], get_data['V1']['eccentricity'][:] ])
template_eccentricity = np.nanmax(eccentricity_vstack, axis=0)

spin1z_vstack = np.vstack([get_data['H1']['spin1z'][:], get_data['L1']['spin1z'][:], get_data['V1']['spin1z'][:]])
template_spin1z = np.nanmax(spin1z_vstack, axis=0)

spin2z_vstack = np.vstack([get_data['H1']['spin2z'][:], get_data['L1']['spin2z'][:], get_data['V1']['spin2z'][:]])
template_spin2z = np.nanmax(spin2z_vstack, axis=0)

template_chirpmass = pycbc.conversions.mchirp_from_mass1_mass2(template_mass1, template_mass2)
template_ecc_chirpmass = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(template_mass1, template_mass2, template_eccentricity)
template_chieff = pycbc.conversions.chi_eff(template_mass1, template_mass2, template_spin1z, template_spin2z)

inj_mass1 = get_data['inj']['mass1'][:]
inj_mass2 = get_data['inj']['mass2'][:]
inj_eccentricity = get_data['inj']['eccentricity'][:]
inj_spin1z = get_data['inj']['spin1z'][:]
inj_spin2z = get_data['inj']['spin2z'][:]

inj_chirpmass = pycbc.conversions.mchirp_from_mass1_mass2(inj_mass1, inj_mass2)
inj_ecc_chirpmass = pycbc.conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(inj_mass1, inj_mass2, inj_eccentricity)
inj_chieff = pycbc.conversions.chi_eff(inj_mass1, inj_mass2, inj_spin1z, inj_spin2z)

diff_dict = {'mass1':template_mass1-inj_mass1,
            'mass2':template_mass2-inj_mass2,
            'eccentricity':template_eccentricity-inj_eccentricity,
            'spin1z':template_spin1z-inj_spin1z,
            'spin2z':template_spin2z-inj_spin2z,
            'chieff':template_chieff-inj_chieff,
            'mchirp':template_chirpmass-inj_chirpmass,
            'mchirp_ecc':template_ecc_chirpmass-inj_ecc_chirpmass
        }

diff_dict.pop('mass1'); diff_dict.pop('mass2'); diff_dict.pop('spin1z'); diff_dict.pop('spin2z')

true_diff_dict = { 'mchirp': search_mchirp - mps_mchirp,
                   'eccentricity': search_ecc - mps_ecc,
                   'chieff': search_chieff - mps_chi_eff,
                   'mchirp_ecc': search_mchirp_ecc - mps_mchirp_ecc
        }


combined_samples = np.array([diff_dict['mchirp'], diff_dict['mchirp_ecc'], diff_dict['eccentricity'], diff_dict['chieff']])
true_values = np.array([true_diff_dict['mchirp'], true_diff_dict['mchirp_ecc'], true_diff_dict['eccentricity'], true_diff_dict['chieff']])
labels = [r"$\Delta\,\mathcal{M}_c^{}$", r"$\Delta\,\mathcal{M}_c^{\rm ecc}$", r"$\Delta\,e_{20}^{}$", r"$\Delta\,\chi^{}_{\rm eff}$"]

quantiles = [
    0.158655, 0.841345,      # 1σ
    0.001350, 0.998650#,      # 3σ
    #0.000000287, 0.999999713 # 5σ
]

levels = [0.5, 0.9]#[0.6827, 0.9973, 0.9999994]

fig = plt.figure(figsize=(7, 5.5))

fig = corner.corner(combined_samples.T, labels = labels, #truths = true_values, 
                bins=60, 
                smooth=1.8, smooth1d=1.8, fill_contours=False, 
                plot_density=True, plot_datapoints = True,
                #quantiles=quantiles, 
                levels=levels,
                fig = fig,  no_fill_contours=True)
axes = np.array(fig.axes).reshape((4,4))

for i in range(4):
    for j in range(4):
        axes[i,j].tick_params(which="both", direction='in')
        #axes[i,j].tick_params(axis='y', direction='in')
        axes[i,j].tick_params(axis='x', labelsize = 16)
        axes[i,j].tick_params(axis='y', labelsize = 16)

color_corner = 'mediumvioletred'
linestyle = 'solid'
lw = 1.5

axes[0,0].axvline(true_diff_dict['mchirp'], color=color_corner, ls = linestyle, lw=lw)

std_mchirp = np.std( diff_dict['mchirp'], ddof=1)
lower_mchirp = np.median( diff_dict['mchirp']) - std_mchirp
upper_mchirp = np.median( diff_dict['mchirp']) + std_mchirp

axes[0,0].axvline( lower_mchirp, color='black', ls='dashed' )
axes[0,0].axvline( upper_mchirp, color='black', ls='dashed' )

axes[1,0].axvline(true_diff_dict['mchirp'], color=color_corner, ls = linestyle, lw=lw)
axes[2,0].axvline(true_diff_dict['mchirp'], color=color_corner, ls = linestyle, lw=lw)
axes[3,0].axvline(true_diff_dict['mchirp'], color=color_corner, ls = linestyle, lw=lw)


axes[1,1].axvline(true_diff_dict['mchirp_ecc'], color=color_corner, ls = linestyle, lw=lw)

std_mchirp_ecc = np.std( diff_dict['mchirp_ecc'], ddof=1)
lower_mchirp_ecc = np.median( diff_dict['mchirp_ecc']) - std_mchirp_ecc
upper_mchirp_ecc = np.median( diff_dict['mchirp_ecc']) + std_mchirp_ecc

axes[1,1].axvline( lower_mchirp_ecc, color='black', ls='dashed' )
axes[1,1].axvline( upper_mchirp_ecc, color='black', ls='dashed' )


axes[2,1].axvline(true_diff_dict['mchirp_ecc'], color=color_corner, ls = linestyle, lw=lw)
axes[3,1].axvline(true_diff_dict['mchirp_ecc'], color=color_corner, ls = linestyle, lw=lw)

axes[2,2].axvline(true_diff_dict['eccentricity'], color=color_corner, ls = linestyle, lw=lw)
axes[3,2].axvline(true_diff_dict['eccentricity'], color=color_corner, ls = linestyle, lw=lw)

axes[3,3].axvline(true_diff_dict['chieff'], color=color_corner, ls = linestyle, lw=lw)

axes[1,0].axhline(true_diff_dict['mchirp_ecc'], color=color_corner, ls = linestyle, lw=lw)

axes[2,0].axhline(true_diff_dict['eccentricity'], color=color_corner, ls = linestyle, lw=lw)
axes[2,1].axhline(true_diff_dict['eccentricity'], color=color_corner, ls = linestyle, lw=lw)
axes[2,2].axhline(true_diff_dict['eccentricity'], color=color_corner, ls = linestyle, lw=lw)

std_eccentricity = np.std( diff_dict['eccentricity'], ddof=1)
lower_eccentricity = np.median( diff_dict['eccentricity']) - std_eccentricity
upper_eccentricity = np.median( diff_dict['eccentricity']) + std_eccentricity

axes[2,2].axvline( lower_eccentricity, color='black', ls='dashed' )
axes[2,2].axvline( upper_eccentricity, color='black', ls='dashed' )

std_chieff = np.std( diff_dict['chieff'], ddof=1)
lower_chieff = np.median( diff_dict['chieff']) - std_chieff
upper_chieff = np.median( diff_dict['chieff']) + std_chieff

axes[3,3].axvline(lower_chieff, color='black', ls='dashed')
axes[3,3].axvline(upper_chieff, color='black', ls='dashed')

axes[3,0].axhline(true_diff_dict['chieff'], color=color_corner, ls = linestyle, lw=lw)
axes[3,1].axhline(true_diff_dict['chieff'], color=color_corner, ls = linestyle, lw=lw)
axes[3,2].axhline(true_diff_dict['chieff'], color=color_corner, ls = linestyle, lw=lw)

fig.subplots_adjust(hspace=0.03, wspace=0.03, left=0.01, bottom=0.01, top=0.99, right=0.99)

fig.savefig("fig6.png")
fig.savefig("fig6.pdf")
