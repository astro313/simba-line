
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
# Specify renderer
# matplotlib.use('Agg')
# Boiler-plate settings for producing pub-quality figures
# 1 point = 1/72 inch
matplotlib.rcParams.update({'figure.figsize'  : (8,5)    # inches
                            ,'font.size'       : 22      # points
                            ,'legend.fontsize' : 16      # points
                            ,'lines.linewidth' : 2       # points
                            ,'axes.linewidth'  : 2.2       # points
                            ,'text.usetex'     : True    # Use LaTeX to layout text
                            ,'font.family'     : "serif" # Use serifed fonts
                            , 'xtick.major.size': 13     # length, points
                            , 'xtick.major.width': 1     # points
                            , 'xtick.minor.size': 8     # length, points
                            , 'xtick.minor.width': 1     # points
                            , 'ytick.major.size': 13     # length, points
                            , 'ytick.major.width': 1     # points
                            , 'ytick.minor.size': 8     # length, points
                            , 'xtick.labelsize': 16
                            , 'ytick.labelsize': 16
                            , 'ytick.minor.width': 1     # points
                            , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"),
                            'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"),
                            'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                            })


import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import caesar
import function as fu
import OBSSMF as obs

nrowmax = 3
fill_between = True

def massFunc(objs, labels, ax, jwind, fill_between=False, decomSFQ=False):
    """
    fill_between: bool
        whether to use fill_between or error bars

    decomSFQ: bool
        whether to plot SF galaxies and Quenched galaxies (based on sSFR) separately.

    """
    for j in range(0, len(objs)):
        for curType in TYPES:
            galpos = np.array([g.pos.d for g in objs[j].galaxies])
            if curType == 'GSMF':
                mass = np.array(
                    [g.masses['stellar'].d for g in objs[j].galaxies])
            elif curType == 'HI':
                mass = np.array([g.masses['HI'].d for g in objs[j].galaxies])
            elif curType == 'H2':
                mass = np.array([g.masses['H2'].d for g in objs[j].galaxies])
            elif curType == 'SFR':
                mass = np.array([g.sfr.d for g in objs[j].galaxies])
            elif curType == 'Halo':
                mass = np.array([h.masses['virial'].d for h in objs[j].halos])
                galpos = np.array([h.pos.d for h in objs[j].halos])
            sfr = np.array([g.sfr.d for g in objs[j].galaxies])
            cent = np.array([g.central for g in objs[j].galaxies])
            volume = objs[j].simulation.boxsize.to('Mpccm').d**3
            redshift = objs[j].simulation.redshift
            ssfr = np.log10(1.e9 * sfr / mass + 10**(-2.7 + 0.3 * redshift))
            ssfrlim = -1.8 + 0.3 * redshift
            ncol = int((len(objs) - 1) / nrowmax + 1)
            icol = int(j / nrowmax)
            irow = int(j % nrowmax)
            if ncol == 1:
                try:
                    ax0 = ax[irow]
                except:
                    ax0 = ax
            else:
                try:
                    ax0 = ax[irow][icol]
                except:
                    ax0 = ax

            x, y, sig = fu.cosmic_variance(
                mass, galpos, objs[j].simulation.boxsize, volume, nbin=18, minmass=7)

            if fill_between:
                ax0.plot(np.log10(x) + j * 0.001, np.log10(y), '--',
                     color=colors[j])
                ax0.fill_between(np.log10(x) + j * 0.001,
                                 np.log10(y) - sig,
                                np.log10(y) + sig,
                                facecolor=colors[j], alpha=0.5,
                                interpolate=True, label=labels[j])
            else:
                ax0.plot(np.log10(x) + j * 0.001, np.log10(y), '--',
                         color=colors[j], label=labels[j])

                ax0.errorbar(np.log10(x)+j*0.001,
                             np.log10(y),
                             yerr=[sig,sig],
                             color=colors[j])

            if decomSFQ:
                # plot SF and Q separately
                # select SF galaxies
                select = [(ssfr > ssfrlim) & (cent == 1)]
                x, y, sig = fu.cosmic_variance(mass[select], galpos[select], objs[
                                               j].simulation.boxsize, volume, nbin=20, minmass=7)
                ax0.plot(np.log10(x) + j * 0.001, np.log10(y), '--',
                         color='b', label=labels[j] + '-SF')

                # Select Quenched galaxies
                select = [(ssfr <= ssfrlim) & (cent == 1)]
                x, y, sig = fu.cosmic_variance(mass[select], galpos[select], objs[
                                               j].simulation.boxsize, volume, nbin=20, minmass=7)
                ax0.plot(np.log10(x) + j * 0.001, np.log10(y), '--',
                         color='r', label=labels[j] + '-Q')

    smf = obs.SMF(redshift)
    if smf.cond:
        ax0.errorbar(smf.x, smf.y, yerr=smf.yerr, fmt='x',
                     label=smf.name, zorder=100, color='k')
    ax0.legend(loc='lower left', fontsize=16)
    # ax0.annotate('z=%g' % np.round(objs[j].simulation.redshift, 1), xy=(
    #     0.8, 0.75), xycoords='axes fraction', size=12, bbox=dict(boxstyle="round", fc="w"))
    plt.title(r'$z$ = %g' % np.round(objs[j].simulation.redshift, 1))
    return None


colors = ('b', 'g', 'm', 'c', 'k')
# simName = ['m25n256', 'm25n1024', 'm50n512', 'm50n1024', 'm100n1024']
simName = ['m25n1024', 'm50n1024', 'm100n1024']
labels = ['Simba-25', 'Simba-50', 'Simba-100']

sims = []
for ii, nn in enumerate(simName):
    if nn == 'm25n1024':
        caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/' + nn + '/s50_new/Groups/' + nn + '_036.hdf5'

    else:
        caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/' + nn + '/s50/Groups/' + nn + '_036.hdf5'
    sims.append(caesar.load(caesarfile, LoadHalo=False))

ncol = 1
nrow = 1

fig, ax = plt.subplots()
TYPES = ['GSMF']
massFunc(sims, labels, ax, ii, fill_between=fill_between, decomSFQ=False)


plt.minorticks_on()
plt.xlim(7,11.5)
plt.xlabel(r'$\log M_* [M_\odot]$',fontsize=18)
plt.ylabel(r'$\log \Phi$ [Mpc$^{-3}$]',fontsize=18)
plt.subplots_adjust(hspace=.0)

if fill_between:
    figname = 'gsmf_test_fill.pdf'
else:
    figname = 'gsmf_test.pdf'
plt.savefig(figname, bbox_inches='tight')
