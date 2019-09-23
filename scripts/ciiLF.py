
"""

LF contains info. about:
* how the population of galaxies change as a function of cosmic time.
* the volume density of these sources as a function of luminosity

* primordial density fluctuations
* processes that destroy/create galaxies
* processes that change one type of galaxy into another (e.g., merger, stripping)
* processes that transfer mass into light


"""

def setup():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns
    # Specify renderer
    # matplotlib.use('Agg')

    # Boiler-plate settings for producing pub-quality figures
    # 1 point = 1/72 inch
    from cycler import cycler
    matplotlib.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

    matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                                , 'font.size': 22      # points
                                , 'legend.fontsize': 16      # points
                                , 'lines.linewidth': 1.5       # points
                                , 'axes.linewidth': 1.5       # points
                                , 'text.usetex': True    # Use LaTeX to layout text
                                , 'font.family': "serif"  # Use serifed fonts
                                , 'xtick.major.size': 10     # length, points
                                , 'xtick.major.width': 1.5     # points
                                , 'xtick.minor.size': 6     # length, points
                                , 'xtick.minor.width': 1     # points
                                , 'ytick.major.size': 10     # length, points
                                , 'ytick.major.width': 1.5     # points
                                , 'ytick.minor.size': 6     # length, points
                                , "xtick.minor.visible": True
                                , "ytick.minor.visible": True
                                , 'font.weight': 'bold'
                                , 'ytick.minor.width': 1     # points
                                , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                                })


    import matplotlib.pyplot as plt
    import pandas as pd
    import sys
    import os


def get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m50n1024/s50/Groups/m50n1024_036.hdf5'):
    """
    Ran on Rusty, where the caesar file is, and caesar is installed
    """
    import caesar      # module purge; module load gcc python3
    bubu = caesar.load(caesarfile, LoadHalo=False)
    volume = bubu.simulation.boxsize.to('Mpccm').d**3
    print(volume)
    return volume


def mass_function(mass, vol_Mpc, nbin=30, minmass=6):
    """

    compute LF or MF

    Parameters
    ----------

    minmass: float or int
        smallest bin = 10^minmass

    """
    maxmass = np.log10(np.nanmax(mass))
    bin_    = np.logspace(minmass, maxmass, nbin+1)
    dM      = ( maxmass - minmass ) / nbin
    logx    = (np.log10(bin_) + dM/2.)[:nbin]
    hist    = np.histogram( mass, bins=bin_, range=(bin_.min(),bin_.max()) )[0]
    y       = hist / vol_Mpc / dM         # Normalize to volume and bin size
    return logx, y, bin_, dM

def schechter_fit(logL, logphi):

    """


    Returns
    -------
    alpha: float
        faint-end slope

    logLstar: float
        characteristic luminosity at noralization point [Lsun]

    logphistar: float
        normalization [Mpc^-3/h^3]

    """
    return alpha, logLstar, logphistar

def schetcher(phistar, Lstar, alpha, Lmin=5., Lmax=9., nbin=100):
    """
    Eqn 20 of Popping+16

    NOTE
    ----
    see also
        https://pythonhosted.org/Astropysics/coremods/models.html
        http://www.simondriver.org/Teaching/AS3011/AS3011_11.pdf

    """
    L = 10**np.linspace(Lmin, Lmax, nbin)
    phi =  phistar * (L/Lstar)**(alpha + 1) * np.exp(-L/Lstar)
    return L, phi


def oplot_popping16_Fig9(fig, ax, nbin=30, alpha=-1.77, logLstar=7.80, logphistar=-2.95):
    # for z=6
    x, y = schetcher(10**logphistar, 10**logLstar, alpha, nbin=nbin)
    ax.plot(np.log10(x), np.log10(y), '--', color='darkorange', label='Popping+16')
    return fig, ax

def oplot_popping19(fig, ax, fname='../../literature/popping19_ciiLF_z6/z6.txt'):
    ### Lsun, /cMpc^3
    x, y = np.loadtxt(fname, unpack=True)
    ax.plot(x, np.log10(y), '--', color='indigo', label='Popping+19')
    return fig, ax


def plot_this_work(mass, cvolume, fig=None, ax=None, nbin=10,
                  errorbar=False, method='jackknife', fill_between=True,
                  galpos=None, color='b', label='This Work',
                  minmass=6.2, saveFig=True):

    import matplotlib.pyplot as plt

    if fig is None and ax is None:
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        fig, ax = oplot_popping16_Fig9(fig, ax)
        fig, ax = oplot_popping19(fig, ax, fname='z6.txt')
        fig, ax = plot_Hemmati17_Capak15(fig, ax, fname='Hemmati17_z6.txt')

    if errorbar:
        if method == 'jackknife':
            assert galpos is not None
            import function as fu
            import OBSSMF as obs
            x, y, sig = fu.cosmic_variance(mass, galpos,
                                          cvolume**(1/3.),
                                          cvolume,
                                          nbin=nbin, minmass=minmass)
            if fill_between:
                ax.plot(np.log10(x), np.log10(y),
                        color=color)

                ax.fill_between(np.log10(x),
                                np.log10(y) - sig,
                                np.log10(y) + sig,
                                facecolor=color, alpha=0.3,
                                interpolate=True, label=label)
            else:
                ax.plot(np.log10(x), np.log10(y),
                         color=color, label=label)

                ax.errorbar(np.log10(x),
                             np.log10(y),
                             yerr=[sig,sig],
                             color=color)
        else:
            raise NotImplementedError

    else:
        logx, y, _, _ = mass_function(mass, cvolume, nbin=nbin, minmass=minmass) # np.log10(mass.min()))
        ax.plot(logx, np.log10(y), label='This Work')


    if saveFig:
        ax.set_ylabel(r'$\log \Phi$ [h$^3$ Mpc$^{-3}$]',fontsize=16)
        ax.set_xlabel(r'$\log L_{\rm [CII]} [L_\odot]$',fontsize=16)
        ax.set_ylim([-6.3, 0.2])
        # ax.set_xlim([5.5, 9.])
        plt.legend(loc='best', fontsize=13)
        plt.title('[CII] LF')
        plt.minorticks_on()
        plt.subplots_adjust(left=0.1, hspace=.0)
        plt.savefig('./CIILF.pdf',bbox_inches='tight')

    return fig, ax


def plot_Hemmati17_Capak15(fig, ax, fname='../../literature/Hemmati17_z6.txt'):
    ### Lsun, /cMpc^3,
    logLcii, phi, phiup, philow = np.loadtxt(fname, unpack=True)
    ax.errorbar(logLcii, np.log10(phi),
                yerr=[0.424 * (phi - philow)/phi, 0.424 * (phiup - phi)/phi],
                label='Hemmati+17',
                alpha=0.95,
                marker='o',
                color='green',
                markeredgecolor='black',
                markeredgewidth=1.0,
                linestyle='None'
                )
    return fig, ax


def M_weighted_pos(pos, mgas):
    weightedpos = np.sum((mgas) * pos) / np.sum((mgas))
    return weightedpos


def get_xyz_from_sigamePD(galnameList, path='/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z6_data_files/particle_data/sim_data_m2550100_7675/', physical=True, redshift=5.93):

    import pandas as pd
    import glob
    import os.path

    pdlist = glob.glob(path+'*.gas')
    pdlist = [os.path.basename(i) for i in pdlist]
    starti = pdlist[0].find('_') + 1
    prefix = pdlist[0][:starti]
    endi = [i.find('_sim') for i in pdlist]

    galn = [n[starti: endi[i]] for i, n in enumerate(pdlist)]

    galpos = []
    missing = []
    for ggg in galnameList:
        if ggg in galn:
            gasProp = pd.read_pickle(path + prefix + ggg + '_sim.gas')
            # starProp = pd.read_pickle(path + prefix + ggg + '_sim.star')
            if physical:                 # physical kpc
                x = M_weighted_pos(gasProp['x'].values, gasProp['m'])
                y = M_weighted_pos(gasProp['y'].values, gasProp['m'])
                z = M_weighted_pos(gasProp['z'].values, gasProp['m'])

                gp = [x, y, z]
            else:             # ckpc
                x = M_weighted_pos(gasProp['x'].values*(1+redshift), \
                                           gasProp['m'])
                y = M_weighted_pos(gasProp['y'].values*(1+redshift), \
                                           gasProp['m'])
                z = M_weighted_pos(gasProp['z'].values*(1+redshift), \
                                           gasProp['m'])
                gp = [x, y, z]

            galpos.append(gp)

        else:
            missing.append(ggg)


    if len(missing) > 0:
        import pdb; pdb.set_trace()

    assert len(galpos) == len(galnameList)
    return np.array(galpos)



if __name__ == '__main__':

    setup()
    size = '2550100'


    if size == '25':
        cvolume = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50_new/Groups/m25n1024_036.hdf5')
    elif size == '50':
        cvolume = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m50n1024/s50/Groups/m50n1024_036.hdf5')
    elif size == '100':
        cvolume = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m100n1024/s50/Groups/m100n1024_036.hdf5')
    elif size == '2550100':
        try:
            cvolume25 = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50_new/Groups/m25n1024_036.hdf5')
            cvolume50 = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m50n1024/s50/Groups/m50n1024_036.hdf5')
            cvolume100 = get_vol(caesarfile = '/mnt/ceph/users/daisyleung/simba/sim/m100n1024/s50/Groups/m100n1024_036.hdf5')
        except ImportError:
            cvolume25 = 49692.77936087925
            cvolume50 = 397542.234887034
            cvolume100 = 3180337.879096272

    else:
        raise NotImplementedError("....")

    import pandas as pd
    import numpy as np
    galgal = pd.read_pickle('combinedPD.pkl')
    mass = galgal['L_CII']

    if size in ['25', '50', '100']:
        plot_this_work(mass, cvolume, nbin=10, minmass=6.2)

    elif size == '2550100':
        #e.g., for m2550100...
        # 1867 from m25
        # 4350 from m50
        # 1458 from m100

        # ad-hoc way for now
        mass25 = mass[:1868]
        mass50 = mass[1868: 1868+4350+1]
        mass100 = mass[1868+4350+1: ]

        # -- no error bars ---
        # fig, ax = plot_this_work(mass25, cvolume25, nbin=10, minmass=4.5, saveFig=False, label='This Work; m25n1024', color='r')
        # fig, ax = plot_this_work(mass50, cvolume50, fig=fig, ax=ax, nbin=10, minmass=6.2, saveFig=False, label='This Work; m50n1024', color='b')
        # plot_this_work(mass100, cvolume100, fig=fig, ax=ax, nbin=10, minmass=7., label='This Work; m100n1024', color='k')

        # add Jackknife error bars; need to run on RUSTY!!!
        galpos = get_xyz_from_sigamePD(galgal['galnames'].values,
                                       path='/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z6_data_files/particle_data/sim_data_m2550100_7675/',
                                       physical=True,   # doesn't actually matter
                                       redshift=5.93)

        galpos /= 1.e3          # cMpc, doesn't actually matter, but just to be consistent w/ other codes

        # ad-hoc
        galpos25 = galpos[:1868]
        galpos50 = galpos[1868: 1868+4350+1]
        galpos100 = galpos[1868+4350+1: ]

        fig, ax = plot_this_work(mass25, cvolume25, nbin=10,
                                 errorbar=True, galpos=galpos25,
                                 minmass=4.5, saveFig=False,
                                 label='Simba-25', color='r')

        fig, ax = plot_this_work(mass50, cvolume50, fig=fig, ax=ax,
                                 nbin=10, errorbar=True, galpos=galpos50,
                                 minmass=6.2, saveFig=False,
                                 label='Simba-50', color='b')

        plot_this_work(mass100, cvolume100, fig=fig, ax=ax, nbin=10,
                       errorbar=True, galpos=galpos100, minmass=7.,
                       label='Simba-100', color='k')

    else:
        raise NotImplementedError("....")
