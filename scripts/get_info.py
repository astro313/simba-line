import caesar
import numpy as np
import cPickle as pickle
import socket
import os
import sys
from readgadget import readsnap

def info(obj, snapFile, top=None, savetxt=False):
    """
    Customized version to print more info on the CAESAR extracted galaxies.

    Parameters
    ----------
    obj : :class:`main.CAESAR`
        Main CAESAR object.
    top : int
        Number of objects to print.

    """

    group_list = obj.galaxies
    nobjs = len(group_list)

    if top > nobjs or top is None:
        top = nobjs

    if obj.simulation.cosmological_simulation:
        time = 'z=%0.3f' % obj.simulation.redshift
    else:
        time = 't=%0.3f' % obj.simulation.time

    output  = '##\n'
    output += '## Largest %d Galaxies \n' % (top)
    if hasattr(obj, 'data_file'): output += '## from: %s\n' % obj.data_file
    output += '## %d @ %s' % (nobjs, time)
    output += '##\n'
    output += '##\n'

    cnt = 1

    output += '## ID      Mstar     Mgas      MBH    fedd    SFR [Msun/yr]      SFRSD [Msun/yr/kpc^2]    SFRSD_r_stellar_half_mass [Msun/yr/kpc^2]    gasSD [Msun/pc^2]    r_baryon   r_gas      r_gas_half_mass      r_stellar    r_stellar_half_mass    Z_sfrWeighted [/Zsun]    Z_massWeighted [/Zsun]     Z_stellar [/Zsun]     T_gas_massWeighted    T_gas_SFRWeighted   fgas   nrho      Central\t|  Mhalo     HID\n'
    output += '## ----------------------------------------------------------------------------------------\n'

    h = obj.simulation.hubble_constant
    bhmdot = readsnap(snapFile, 'BH_Mdot','bndry',suppress=1)*1.e10/h/3.08568e+16*3.155e7 # in Mo/yr

    # debug mbh:
    # mmm = [g.masses['bh'] for g in obj.galaxies]
    # mmm[:10]

    for ii, o in enumerate(group_list):
        phm, phid = -1, -1

        if o.halo is not None:
            phm, phid = o.halo.masses['total'], o.halo.GroupID

        try:
            bhmdots = [bhmdot[k] for k in o.bhlist]
            bm = o.masses['bh']
            imax = np.argmax(bm)
            try:
                bm = bm[imax]
                bmdot = bhmdots[imax]        # only the massive BH particle matters.
            except:
                bm = bm
                bmdot = bhmdots
            frad = 0.1
            mdot_edd = 4*np.pi*6.67e-8*1.673e-24/(frad*3.e10*6.65245e-25) * bm * 3.155e7 # in Mo/yr
            fedd = bmdot / mdot_edd
            fedd = fedd[0].value
        except:
            bm = 0
            fedd = 0

        # print bm, fedd
        # import pdb; pdb.set_trace()

        output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e   %0.3f   %0.3f  %0.3f  %0.2e  %0.2e  %0.3f  %0.2e  %s\t|  %0.2e  %d \n' % \
                  (o.GroupID, o.masses['stellar'], o.masses['gas'],
                   bm,
                   fedd,
                   o.sfr,
                   o.sfr/np.pi/o.radii['gas']**2,
                   o.sfr/np.pi/o.radii['stellar_half_mass']**2,
                   o.masses['gas']/np.pi/(o.radii['gas']*1.e3)**2,
                   o.radii['baryon'],
                   o.radii['gas'],
                   o.radii['gas_half_mass'],
                   o.radii['stellar'],
                   o.radii['stellar_half_mass'],
                   o.metallicities['sfr_weighted']/0.0134,
                   o.metallicities['mass_weighted']/0.0134,
                   o.metallicities['stellar']/0.0134,
                   o.temperatures['mass_weighted'],
                   o.temperatures['sfr_weighted'],
                   o.gas_fraction,
                   o.local_number_density, o.central,
                   phm, phid)

        cnt += 1
        if cnt > top: break

    print(output)

    if savetxt:

        outputFileName = snapFile[50:-5] + '_top' + str(top) + '.txt'
        f = open(outputFileName, 'w')
        f.write(output)
        f.close()

    return output, outputFileName


def setup_cmap(cm='magma'):    # "plasma_r"
    import matplotlib.pyplot as plt
    plt.set_cmap(cm)
    return plt.get_cmap()


def setup_plot():
    import matplotlib
    from cycler import cycler
    # print(matplotlib.matplotlib_fname())

    matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                                , 'font.size': 16      # points
                                , 'legend.fontsize': 10      # points
                                , 'lines.linewidth': 2       # points
                                , 'axes.linewidth': 1       # points
                                , 'axes.prop_cycle': cycler('color', 'bgrcmyk')
                                , 'text.usetex': True
                                , 'font.family': "serif"  # Use serifed fonts
                                , 'xtick.major.size': 13     # length, points
                                , 'xtick.major.width': 1     # points
                                , 'xtick.minor.size': 8     # length, points
                                , 'xtick.minor.width': 1     # points
                                , 'ytick.major.size': 13     # length, points
                                , 'ytick.major.width': 1     # points
                                , 'ytick.minor.size': 8     # length, points
                                , 'xtick.labelsize': 14
                                , 'ytick.labelsize': 14
                                , 'ytick.minor.width': 1     # points
                                , 'font.serif': ("times", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                                })
    cm = setup_cmap()
    return cm


def plot_info(colNumx, colNumy, inFile,
              colNumz=None, zlabel='',
              logx=True, logy=True, logz=True,
              ythreshold=None,
              xthreshold=None,
              xlabel='',
              ylabel='',
              ls='',
              marker='*',
              markersize=10,
              legendFontSize=17,
              cbarLabelSize=17,
              tag='',
              cmap=None,     # 'viridis'
              saveFig=True, savedir='../plots/'):
    """

    colNumx: int
        column to read from the .txt file, used as x-axis in plotting

    colNumy: int
        column to read from the .txt file, used as y-axis in plotting

    colNumz: int
        column to read from the .txt file, used as color in plotting

    inFile: str
        name of input file to read from where we extract global properties of galaxies


    """

    import matplotlib

    colNumx = int(colNumx)
    colNumy = int(colNumy)
    xxx, yyy = np.genfromtxt(inFile, usecols=(colNumx, colNumy), unpack=True)

    if colNumz:
        colNumz = int(colNumz)
        zzz = np.genfromtxt(inFile, usecols=(colNumz))

    import matplotlib.pyplot as plt
    cm = setup_plot()
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)

    print("filter out bad points...")
    if xlabel[0] == 'M':
        bad = xxx == 0.0
        xxx = xxx[~bad]
        yyy = yyy[~bad]
        if colNumz:
            zzz = zzz[~bad]

    if ylabel[0] == 'M':
        bad = yyy == 0.0
        xxx = xxx[~bad]
        yyy = yyy[~bad]
        if colNumz:
            zzz = zzz[~bad]

    if not ythreshold is None:
        # select points with threshold above certain values in yaxis
        good = yyy > ythreshold
        yyy = yyy[good]
        xxx = xxx[good]
        if colNumz:
            zzz = zzz[good]

    if not xthreshold is None:
        good = xxx > xthreshold
        yyy = yyy[good]
        xxx = xxx[good]
        if colNumz:
            zzz = zzz[good]

    if colNumz:
        fig.subplots_adjust(right=0.84, wspace=0.01)

        if logz:
            if zzz.min() == 0.0:
                zzz += 1.e-6

            c = list(np.log10(zzz))

        else:
            c = list(zzz)

        if cmap is None:
            cmap = plt.get_cmap()
        else:
            cmap = matplotlib.cm.get_cmap(cmap)
        normalize = matplotlib.colors.Normalize(vmin=min(c), vmax=max(c))
        colors = [cmap(normalize(value)) for value in c]
        # NUM_COLORS = len(c)
        # ax.set_prop_cycle('color', [cm(1. * i / NUM_COLORS)
        #                             for i in range(NUM_COLORS)])

        ax.scatter(xxx, yyy, color=colors, s=markersize,
                 marker=marker, label='This Work')

    else:
        ax.plot(xxx, yyy, ls=ls, markersize=markersize,
                 marker=marker, color='k',
                 # markeredgecolor='gray',
                 markeredgewidth=0.5, label='This Work')

    if logx:
        # ax.set_xscale('symlog')
        ax.set_xscale('log')

    if logy:
        # ax.set_yscale('symlog')
        ax.set_yscale('log')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.minorticks_on()

    if colNumz:
        cax, _ = matplotlib.colorbar.make_axes(ax)
        cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)

        cbar.ax.tick_params(length=6, labelsize=legendFontSize)
        cbar.set_label(zlabel, fontsize=cbarLabelSize)

    # literature SK
    if xlabel.lower() == "gassd" and "sfrsd" in ylabel.lower():
        litpath = './literature/'

        x, y = [10**0.50, 10**4.0], [10**(-2.85), 10**2.1]
        ax.plot(x, y, linestyle='-', color='b',
                linewidth=2, label="Kennicutt 1998")

        # more from high-z literature
        x0901, y0901 = np.loadtxt(
            litpath + "J0901_KS10_points.txt", unpack=True)  # in log
        x14011, y14011 = np.loadtxt(
            litpath + "J14011_KSpoints2.txt", unpack=True)  # not in log
        xrawle, yrawle = np.loadtxt(
            litpath + "Rawle_KSpoints.txt", unpack=True, usecols=(0, 1))  # in log
        xgn20, xgn20err, ygn20, ygn20err = np.loadtxt(
            litpath + "Hodge_resolvedKS.txt", unpack=True)   # not in log
        xegs, xegserr, yegs, yegserr = np.loadtxt(
            litpath + "Genzel_KSpoints.txt", unpack=True)  # in log

        ax.scatter(10**x0901, 10**y0901, label="J0901 @ z=2.26",
                   color='red', marker='o', s=5, facecolors='none', alpha=0.6)
        ax.scatter(x14011, y14011, label="SMM J14011 @ z=2.56",
                   color='darkblue', marker='^', s=13, facecolors='none', alpha=0.8)
        ax.scatter(10**xrawle, 10**yrawle, label="HLS0918 @ z=5.24",
                   color='purple', marker='v', s=10, facecolors='none', alpha=0.8)
        ax.errorbar(xgn20, ygn20, yerr=ygn20err, xerr=xgn20err,
                    label="GN20 @ z=4.05",
                    color='orange', fmt='s', markersize=4.5,
                    markeredgewidth=0.6, mfc='none', elinewidth=0.5)
        ax.errorbar(10**xegs, 10**yegs, yerr=yegserr, xerr=xegserr,
                    label="EGS13011166 @ z=1.53",
                    color='green', fmt='D', markersize=3.5,
                    markeredgewidth=0.6, mfc='none', zorder=0.5, alpha=0.56, elinewidth=0.5)
        ax.legend(loc='best', ncol=2)

    if not saveFig:
        plt.show(block=False)
    else:
        # fig.subplots_adjust(right=0.84, left=0.1, top=0.85,
        #                     bottom=0.1, hspace=0.2,
        #                     wspace=0.1)
        if not os.path.isdir(savedir):
            os.makedirs(savedir)

        figName = savedir + xlabel + '_' + ylabel
        if colNumz:
            figName += '_' + zlabel
        if len(tag) > 0:
            figName += '_' + tag + '_'
        figName += '.png'
        plt.savefig(figName, bbox_inches='tight')

    return fig, ax


def plot_literature_SK():
    return fig, ax


if __name__ == '__main__':

    if len(sys.argv) > 1:
        debug = sys.argv[1]
        LoadHalo = False
    else:
        debug = False
        LoadHalo = True

    snapRange = [36]    # don't put 036
    zCloudy = 6
    raw_sim_name_prefix = 'snap_m25n1024_'
    name_prefix = 'm25n1024_'

    host = socket.gethostname()
    if 'ursa' in host:
        raw_sim_dir = '/disk01/rad/sim/m25n1024/s50/'
        caesar_dir = '/disk01/rad/sim/m25n1024/s50/Groups/'
        redshiftFile = '/home/rad/gizmo-extra/outputs_boxspace50.info'
        d_data = '/home/dleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'
    elif 'flatironinstitute.org' or 'worker' in host:
        raw_sim_dir = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50/'  # dummy
        caesar_dir = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50/Groups/'
        redshiftFile = '/mnt/ceph/users/daisyleung/simba/gizmo-extra/outputs_boxspace50.info'
        d_data = '/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'

        if 'worker' in host:
            import matplotlib
            matplotlib.use('Agg')

    infile = caesar_dir + name_prefix + '{:0>3}'.format(int(36)) + \
                '.hdf5'
#    obj = caesar.load(infile, LoadHalo=LoadHalo)
    snapFile = raw_sim_dir + raw_sim_name_prefix + '{:0>3}'.format(int(36)) + '.hdf5'

    output, outName = info(obj, snapFile, top=None, savetxt=True)
    savedir='../plots/' + str(outName[:outName.find('.txt')]) + '/'

    fig, ax = plot_info(1, 22+1, inFile=outName, colNumz=19, zlabel='fgas',
                        xlabel='Mstar', ylabel='Mhalo', logz=False,
                        savedir=savedir)

    fig, ax = plot_info(2, 1, inFile=outName, colNumz=3, xlabel='Mgas', \
                        ylabel='Mstar', zlabel='MBH', savedir=savedir)

    fig, ax = plot_info(2, 1, inFile=outName, colNumz=5, xlabel='Mgas', \
                        ylabel='Mstar', zlabel='SFR', logz=False,
                        savedir=savedir)
    # MZR:
    fig, ax = plot_info(1, 16, inFile=outName, colNumz=14, xlabel='Mstar', \
                    ylabel='Zstellar',
                    logy=False,
                    zlabel='Zgas', logz=False,
                    savedir=savedir)

    fig, ax = plot_info(1, 14, inFile=outName, colNumz=19, xlabel='Mstar', \
                    ylabel='Zgas', ythreshold=0.0,
                    logy=False,
                    zlabel='fgas', logz=False, savedir=savedir)

    # FMR: SFR - Z - M*
    fig, ax = plot_info(5, 14, inFile=outName, colNumz=1, xlabel='SFR', \
                    xthreshold=0.1,     # select SFR > 0.1
                    logy=False,
                    ylabel='Zgas', ythreshold=0.0,
                    zlabel='Mstar', savedir=savedir)

    # SFR f_gas
    fig, ax = plot_info(5, 19, inFile=outName, colNumz=1, xlabel='SFR', \
                    xthreshold=0.1,
                    ylabel='fgas', logy=False, zlabel='Mstar', savedir=savedir)

    fig, ax = plot_info(6, 19, inFile=outName, colNumz=9,
                    xlabel='SFRSD [Msun/pc2]',
                    xthreshold=0.0,
                    ylabel='fgas', logy=False, zlabel='Rbaryon [kpc]',
                    logz=False, savedir=savedir)

    # SFRSD - GasSD
    fig, ax = plot_info(8, 6, inFile=outName, colNumz=9,
                        xlabel='gasSD [Msun/pc2]', xthreshold=0.0,
                        ylabel='SFRSD [Msun/yr/kpc2]', ythreshold=0.0,
                        zlabel='Rbaryon [kpc]',
                        logz=False,
                        savedir=savedir)

    # SFRSD_stellar half mass R - GasSD
    fig, ax = plot_info(8, 7, inFile=outName, colNumz=9,
                        xlabel='gasSD [Msun/pc2]', xthreshold=0.0,
                        ylabel='SFRSDrstellarhalfmass [Msun/yr/kpc2]',
                        ythreshold=0.0,
                        zlabel='Rbaryon [kpc]',
                        logz=False,
                        savedir=savedir)

