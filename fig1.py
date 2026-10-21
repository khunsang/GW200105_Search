import matplotlib
matplotlib.use('Agg')
from pycbc import conversions
from pycbc.waveform import bank
import matplotlib.pyplot as plt
import mpl_scatter_density # adds projection='scatter_density'
from cycler import cycler

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

def using_mpl_scatter_density(fig, x, y, colormap):
    #norm = ImageNormalize(vmin=0., vmax=1000, stretch=LogStretch())
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    density = ax.scatter_density(x, y, cmap=colormap)#, norm=norm)
    fig.colorbar(density, label=r'\# of templates')


colormap = 'Spectral'

# Eccentric template bank
fname = "Eccentric_template_bank_GW200105_search.hdf"

tmplt_bank = bank.TemplateBank(fname)

templates = tmplt_bank.table

mass1 = templates.mass1
mass2 = templates.mass2
spin1z = templates.spin1z
spin2z = templates.spin2z
eccentricity = templates.eccentricity
mchirp = conversions.mchirp_from_mass1_mass2(mass1, mass2)
mchirp_eccentric = conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(mass1, mass2, eccentricity)
chi_eff = conversions.chi_eff(mass1, mass2, spin1z, spin2z)


# Paremeters of the best matched template with the event GW200105, see Table II of Phukon et al., arxiv:2512.10803
search_m1 = 9.18
search_m2 = 2.06
search_s1z = -0.08
search_s2z = 0.05
search_ecc = 0.105
search_mchirp = conversions.mchirp_from_mass1_mass2(search_m1, search_m2)
search_mch_ecc = conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(search_m1, search_m2, search_ecc)
search_chi_eff = conversions.chi_eff(search_m1, search_m2, search_s1z, search_s2z)



# Paremeters of the median posterior value for GW200105 from G. Morras et al., Astrophys. J. Lett. 1000, L2 (2026)
mps_m1_source = 11.5
mps_m2_source = 1.5
mps_mchirp_source = 3.33
mps_ecc = 0.145
mps_z = 0.06
mps_chi_eff = 0.004



mps_mch_ecc = conversions.mchirp_eccentric_from_mass1_mass2_eccentricity(mps_m1_source, mps_m2_source, mps_ecc)*(1+mps_z)


fig = plt.figure(figsize=(7,5))
using_mpl_scatter_density(fig, mchirp_eccentric, chi_eff, colormap)
plt.plot(search_mch_ecc, search_chi_eff, marker='*', color='white',  markersize=7, label='Triggering template' )
plt.plot(mps_mch_ecc, mps_chi_eff, marker='D', color='white',  markersize=3.5, label='Median value from MPS 2025' )
plt.text(search_mch_ecc, -0.095-0.095*.1, 'Triggering template', horizontalalignment='center', size='small', color='white' )
plt.text(mps_mch_ecc, 0.02+0.02*.1, 'MPS 2025', horizontalalignment='center', size='small', color='white')
plt.xlabel(r"$\mathcal{M}^{\mathrm{ ecc}}_c [\mathrm{M}^{}_\odot]$")
plt.ylabel(r"$\chi_{\mathrm{eff}}^{}$")
fig.tight_layout()
fig.savefig("fig1a.pdf")
fig.savefig("fig1a.png")





fig = plt.figure(figsize=(7,5))
using_mpl_scatter_density(fig, mchirp, eccentricity, colormap)
plt.plot(search_mchirp, search_ecc, marker='*', color='black',  markersize=7, label='Triggering template' )
plt.plot(mps_mchirp_source*(1+mps_z), mps_ecc, marker='D', color='black',  markersize=3.5, label='Median value from MPS 2025' )
plt.xlabel(r"$\mathcal{M}^{}_c [\mathrm{M}^{}_\odot]$")
plt.ylabel(r"$e_{20}^{}$")
plt.yscale('log')
fig.tight_layout()
fig.savefig("fig1b.pdf")
fig.savefig("fig1b.png")

