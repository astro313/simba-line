
import pylab as plt
import numpy as np
from scipy import stats


def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))


def plotmedian(x, y, yflag=[], c='k', ltype='--', lw=3, stat='median', ax='plt', bins=8, label=None, pos=None, boxsize=-1):
    if len(yflag) != len(x):
        # print 'Plotmedian: No flag provided, using all values'
        xp = x
        yp = y
    else:
        xp = x[yflag]
        yp = y[yflag]
    # bins<0 sets bins such that there are equal numbers per bin
    if bins < 0: bin_edges = histedges_equalN(xp, -bins)
    else: bin_edges = np.arange(0.999 * min(xp), 1.001 * max(xp), (max(xp) - min(xp)) / (bins))
    if bins < 0:    bin_means, bin_edges, binnumber = stats.binned_statistic(
        xp, yp, bins=bin_edges, statistic=stat)
    else: bin_means, bin_edges, binnumber = stats.binned_statistic(xp, yp, bins=bins, statistic=stat)
    bin_cent = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    ax.plot(bin_cent, bin_means, ltype, lw=lw, color=c, label=label)

    if boxsize > 0:
     # determine cosmic variance over 8 octants, plot errorbars
        if len(yflag) != len(x):
            posp = pos

        else:
            posp = pos[yflag]
        pos = np.floor(posp / (0.5 * boxsize)).astype(np.int)
        gal_index = pos[:, 0] + pos[:, 1] * 2 + pos[:, 2] * 4
        bin_oct = np.zeros((abs(bins), 8))
        for i0 in range(8):
            xq = xp[gal_index == i0]
            yq = yp[gal_index == i0]

            bin_oct[:, i0], bin_edges, binnumber = stats.binned_statistic(xq, yq, bins=bin_edges, statistic=stat)
        bin_oct = np.ma.masked_invalid(bin_oct)
        var = np.ma.std(bin_oct, axis=1)
        ax.errorbar(bin_cent, bin_means, yerr=[
                    var, var], fmt='o', linewidth=lw, color=c)
    elif boxsize == -1:
        var = []
        for i0 in range(len(bin_edges) - 1):
            xq = xp[(xp > bin_edges[i0]) & (xp < bin_edges[i0 + 1])]
            yq = yp[(xp > bin_edges[i0]) & (xp < bin_edges[i0 + 1])]
            var.append(np.ma.std(yq - np.mean(yq)))
        # print 'bins',bin_cent
        # print 'variance',var
        ax.errorbar(bin_cent, bin_means, yerr=[
                    var, var], fmt='o', linewidth=lw, color=c)
    else:
        var = np.zeros(len(bin_cent))
    return bin_cent, bin_means, var


def runningmedian(x, y, xlolim=-1.e20, ylolim=-1.e20, bins=10, stat='median'):
    xp = x[(x > xlolim) & (y > ylolim)]
    yp = y[(x > xlolim) & (y > ylolim)]
    if bins < 0:    # bins<0 sets bins such that there are equal numbers per bin
        bin_edges = histedges_equalN(xp,-bins)
        bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bin_edges,statistic=stat)
    else:
        bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bins,statistic=stat)
        bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ymed = []
        ymean = []
        ysigma = []
        for i in range(0,len(bin_edges[:-1])):
                xsub = xp[xp>bin_edges[i]]
                ysub = yp[xp>bin_edges[i]]
                ysub = ysub[xsub<bin_edges[i+1]]
                ymed.append(np.median(10**ysub))
                ymean.append(np.mean(10**ysub))
                ysigma.append(np.std(10**ysub))
        if stat=='median': ymean = np.asarray(ymed)
        else: ymean = np.asarray(ymean)
        ysiglo = np.maximum(ymean-ysigma,ymean*0.1)
        ysiglo = np.log10(ymean)-np.log10(ysiglo)
        ysighi = np.log10(ymean+ysigma)-np.log10(ymean)
        ymean = np.log10(ymean)
        # print bin_cent,ymean,ysiglo,ysighi
        # plt.plot(bin_cent,ymed,'ro',ms=12,color='c')
        # plt.plot(bin_cent,ymean,'--',lw=3,color='m')
        # plt.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro')
    return bin_cent,ymean,ysiglo,ysighi

