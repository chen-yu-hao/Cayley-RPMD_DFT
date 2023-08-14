'''
Author: Wenbin FAN
First Release: Oct. 30, 2019
Modified: 2022-07-11 21:19:44 Wenbin FAN @FDU
Verision: 1.8
[Intro]
Plot
# (1) the overlap of each xi window,
# (2) the variance in each trajectory,
# (3) potential of mean force,
# (4) recrossing factor,
# (5) the evolution of normalized population of xi,
# (6) the force constant,
# (7) the evolution of the potential of mean force (PMF),
# (8) the evolution of free energy and corresponding reaction coordinates,
# (9) the evolution of xi,
# (10) the deviation of xi,
for a single task.
[Usage]
run `python <this file name>` then you will be ask to input a path containing a result for a single task.
Then you'll get all four figures above if your task ended normally.
Attention, please! The former figures will be deleted when the program started running.
[Contact]
Mail: fanwenbin@shu.edu.cn, langzihuigu@qq.com
WeChat: 178-6567-1650
Thanks for using this python program, with the respect to my supervisor Prof. Yongle!
Great thanks to Yang Hui.
[Bug Fixing]
V1.6:
1) compute variance in each traj, pump the large jump.
2) show upper bound of force constants.
3) modified the gradient color of variance.
4) stretch the figure size of overlapping.
V1.5:
1) enhance the stability
2) add more figure
V1.4:
1) plot 10 UI and PMF figures in evolution,
2) plot 3D version of UI and PMF
3) optimize the reading procedure
V1.3:
1) read information from input file,
2) add a plot for force constant.
3) add a plot for potential of mean force.
V1.2:
1) PMF: Plot range modified.
'''
import argparse
parser = argparse.ArgumentParser()
parser.description='please enter two parameters t for temperature n for number of beads...'
parser.add_argument("-t", "--t", help="this is parameter t", dest="T", type=str, default="1000")
parser.add_argument("-n", "--n", help="this is parameter n",dest="N",  type=int, default="64")
parser.add_argument("-i", "--input", help="path of the input.py",dest="I",  type=str, default="RPMDrate/input.py")
parser.add_argument("-R", "--RPMDpath", help="this is parameter n",dest="R",  type=str, default="RPMDrate/")
args = parser.parse_args()
try:
    args.T=int(args.T)
except:
    args.T=float(args.T)
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import numpy as np
import pandas as pd

color = ['#00447c', '#ae0d16', '#47872c', '#800964']
# SHU Blue, Weichang Red, Willow Green, SHU Purple
# This color scheme can be easily obtained on the official website `vi.shu.edu.cn`.
Tcolor1 = [0., 68., 124.]
Tcolor2 = [174., 13., 22.]

def input_path():
    path = input('Please input the folder with `submitting script` and your input file: \n')
    return path


def myEnding():
    print('\n'
          '[INFO] All jobs have been done normally! \n'
          '       Any question please contact `fanwenbin[at]shu.edu.cn`\n')
    return


def clearFolder(path):
    if os.path.exists(os.path.join(figPath, path)):
        for fileName in os.listdir(os.path.join(figPath, path)):
            os.remove(os.path.join(figPath, path, fileName))
    else:
        os.makedirs(os.path.join(figPath, path))


def plot_parameters(title, width=4):
    print('[INFO] Plotting {}! '.format(title))
    plt.figure(figsize=(width, 3))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"  # The font closet to Times New Roman
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.left'] = True
    plt.rcParams['xtick.top'] = True
    plt.minorticks_on()  # Turn on minor ticks

    file = '{}.png'.format(title)
    if os.path.exists(file):
        os.remove(file)


def plot_save(name):
    plt.tight_layout()
    plt.savefig(os.path.join(figPath, '%s_%s_'%(str(args.T),str(args.N))+name + '.png'), format='png', dpi=600)  # .svg is recommended!
    plt.clf()
    plt.close()


def myFormatter():
    sciFormatter = ticker.ScalarFormatter(useMathText=True)
    sciFormatter.set_scientific(True)
    sciFormatter.set_powerlimits((-1, 1))
    return sciFormatter

def my_gaussian(x, xav, xav2):
    y = (1.0 / (np.sqrt(2.0 * np.pi * xav2))) * np.exp(-(x - xav) ** 2 / (2.0 * xav2))
    return y


def plot_overlap():
    title = 'Overlap'
    plot_parameters(title)
    plt.figure(figsize=(9, 3))

    resolution = 2000
    extend = 0.03  # 3E-2

    xiMin = np.min(xi_list)
    xiMax = np.max(xi_list)
    length = len(xi_list)

    x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)
    y_sum = np.zeros((resolution))  # Total density line

    maxPop = 0  # maximum of the summation of all population
    for i in range(length):
        # Gaussian smearing
        xav, xav2 = umbInfo[3:, NtrajEff - 1, i]
        if xav2 < 1E-10:  # xav2 is zero!
            plt.axvline(xi_list[i], ls='--', c='red', lw=0.2)
            print('[ERROR] Variance at `xi = {}` is ZERO! '.format(xi_list[i]))
        else:
            y_new = my_gaussian(x_new, xav, xav2)

            # Find biggest population
            if max(y_new) > maxPop:
                maxPop = max(y_new)

            y_sum += y_new  # sum all population
            if xav2 > 1.0E-4:
                print("[WARNING] May be too various in xi = {}! ".format(xi_list[i]))
                plt.plot(x_new, y_new, lw=1, c=color[1], alpha=0.8)
            else:
                plt.plot(x_new, y_new, lw=0.5, c=color[0], alpha=.3)

    # Plot summation and difference
    plt.plot(x_new, y_sum, lw=1, c=color[0], label=mylabel)  # label='Summation of all populations')  # SHU Blue

    # plt.xlabel('Reaction Coordinate / Å')
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('Population')

    plt.xlim(xiMin - extend, xiMax + extend)
    # plt.ylim(0, maxPop*1.1)
    plt.ylim(0, max(y_sum) * 1.2)

    plt.yticks([])  # No ticks and labels in y axis

    plt.legend(loc='upper left')
    plot_save(title)


def plot_variance():
    if NtrajEff == 1:
        return

    title = 'Variance'
    plot_parameters(title)

    xiMin = np.min(xi_list)
    xiMax = np.max(xi_list)
    length = len(xi_list)

    timeMax = 0.0
    timeMin = 1E5
    jumpCount = 0

    for i in range(length):
        xivar = umbInfo[4, :, i]
        timeEvolution = umbInfo[2, :, i]

        timeStep = delta
        timeEvolution = [x * timeStep for x in timeEvolution]  # 0.1 fs to 1 ns

        if xivar[-1] > 5E-5:
            print('       traj. at xi = {} may be too various! '.format(xi_list[i]))

        # largeJumpJudge = np.var(xivar)*1E13
        # if largeJumpJudge > 10:
        #     print('{:.3f}\t{:.2f}'.format(xi_list[i], largeJumpJudge))

        for j in range(len(xivar)):
            if j > 0:
                varCurrent = (umbInfo[1, j, i] - umbInfo[1, j - 1, i]) / (umbInfo[2, j, i] - umbInfo[2, j - 1, i]) - \
                             ((umbInfo[0, j, i] - umbInfo[0, j - 1, i]) / (
                                         umbInfo[2, j, i] - umbInfo[2, j - 1, i])) ** 2
                if varCurrent / umbInfo[4, j, i] - 1 > 1:
                    if jumpCount == 0:
                        print('[INFO] Large jump: xi, step/10000, Current var, var')
                    print('{:.3f}\t{:>5d}\t{:.3e}\t{:.3e}'.format(xi_list[i], int(umbInfo[2, j, i] / 10000),
                                                                  varCurrent, umbInfo[4, j, i]))
                    jumpCount += 1

        # xivarDelta = []
        # for i in range(len(xivar)-1):
        #     xivarDelta.append(np.abs(xivar[i+1] - xivar[i]))
        # x = range(len(xivar))
        # plt.yscale('log')

        # # Shifted
        # for i in range(len(xivar)):
        #     xivar[i] = xivar[i] - xivar[0]

        # define transition color # from SHU blue (xi = 0) to Weichang red (xi = 1)
        # tscolor = (int(((Tcolor2[0] - Tcolor1[0]) * i / length + 0)) / 255.0,
        #            (int((Tcolor2[1] - Tcolor1[1]) * i / length + 68)) / 255.0,
        #            (int((Tcolor2[2] - Tcolor1[2]) * i / length + 124)) / 255.0)
        tscolor = (int((255. * i / length)) / 255.0,
                   0.,
                   (int(-255. * i / length + 255)) / 255.0)
        # colorList = list(tscolor)
        # colorPrint = []
        # for j in range(3):
        #     colorPrint.append(int(colorList[j] * 255))
        # print('{:.3f}\t{:}'.format(xi_list[i], colorPrint))

        # color = (np.random.rand(), np.random.rand(), np.random.rand())
        plt.plot(timeEvolution, xivar, lw=1, c=tscolor, alpha=0.5)

        if max(timeEvolution) > timeMax:
            timeMax = max(timeEvolution)
        if min(timeEvolution) < timeMin:
            timeMin = min(timeEvolution)

    plt.xlabel('Time (ps)')
    plt.ylabel('Variance')

    # print(umbInfo[2, 0, 0] * delta * 1E-3, umbInfo[2, -1, 0] * delta * 1E-3)
    plt.xlim(timeMin, timeMax)
    # 2 is 3rd column, -1 is the time in last frame, 0 is random (other number is ojbk).

    # Scientific notation for y axis
    # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    plt.gca().yaxis.set_major_formatter(formatter)

    plot_save(title)


def plot_variance_diff():
    if NtrajEff == 1:
        return

    title = 'Variance_diff'
    plot_parameters(title)

    xiMin = np.min(xi_list)
    xiMax = np.max(xi_list)
    length = len(xi_list)

    timeMax = 0.0
    timeMin = 1E5
    jumpCount = 0

    v = []

    for i in range(length):
        av = umbInfo[0, :, i]
        av2 = umbInfo[1, :, i]
        timeEvolution = umbInfo[2, :, i]

        Nsample = len(av)
        variance = np.zeros(Nsample)
        for j in range(Nsample):
            if j == 0:
                dmean = av[0] / timeEvolution[0]
                variance[0] = av2[0] / timeEvolution[0] - dmean * dmean
            else:
                dt = timeEvolution[j] - timeEvolution[j-1]
                dmean = (av[j] - av[j-1])  / dt
                variance[j] = (av2[j] - av2[j-1]) / dt - dmean * dmean
                # print(i, xi_list[i], variance[j], dt, dmean, av[j], av2[j])

        v.append(variance)


        tscolor = (int((255. * i / length)) / 255.0,
                   0.,
                   (int(-255. * i / length + 255)) / 255.0)

        timeStep = delta
        timeEvolution = [x * timeStep for x in timeEvolution]  # 0.1 fs to 1 ns

        plt.plot(timeEvolution, variance, lw=0.2, c=tscolor, alpha=0.4)

        if max(timeEvolution) > timeMax:
            timeMax = max(timeEvolution)
        if min(timeEvolution) < timeMin:
            timeMin = min(timeEvolution)

    plt.xlabel('Time (ps)')
    plt.ylabel('Variance')

    # print(umbInfo[2, 0, 0] * delta * 1E-3, umbInfo[2, -1, 0] * delta * 1E-3)
    plt.xlim(timeMin, timeMax)
    # 2 is 3rd column, -1 is the time in last frame, 0 is random (other number is ojbk).

    # Scientific notation for y axis
    # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    plt.gca().yaxis.set_major_formatter(formatter)

    plot_save(title)

    title = 'Variance_diff_box'
    plot_parameters(title, width=8)
    plt.boxplot(v, positions=xi_list, widths=0.003, labels=None, sym='x')
    plt.xlim(xiMin, xiMax)

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    plt.gca().yaxis.set_major_formatter(formatter)

    plot_save(title)

    return


def plot_var_evolution():
    if NtrajEff == 1:
        return

    title = 'Variance in each windows'
    clearFolder('Variances')
    plot_parameters(title)
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))

    xiMin = np.min(xi_list)
    xiMax = np.max(xi_list)
    length = len(xi_list)

    timeMax = np.max(umbInfo[2, :, :]) * 1E-3 * delta

    for i in range(length):

        fig = plt.figure()
        gs = gridspec.GridSpec(2, 2)

        me = fig.add_subplot(gs[0, 0])  # mean evolution
        ms = fig.add_subplot(gs[0, 1])  # mean scatter
        ee = fig.add_subplot(gs[1, 0])  # each evolution
        es = fig.add_subplot(gs[1, 1])  # each scatter

        xiMean = np.add(np.divide(umbInfo[0, :, i], umbInfo[2, :, i]), -xi_list[i])
        varMean = umbInfo[4, :, i]

        xiSum = umbInfo[0, :, i]
        varSum = umbInfo[1, :, i]

        stepCount = umbInfo[2, :, i]
        timeCount = np.multiply(stepCount, delta)
        timeCount_eff = 0
        for k, kkk in enumerate(timeCount[1:]):
            if timeCount[k + 1] > timeCount[k]:
                timeCount_eff += 1
            else:
                break

        xiSep = np.zeros(Ntraj)
        varSep = np.zeros(Ntraj)

        for j in range(Ntraj):
            if j == 0:
                xiSep[j] = umbInfo[3, 0, i]
                varSep[j] = varMean[j]
            elif j > 0:
                stepsCurrent = stepCount[j] - stepCount[j - 1]
                varSep[j] = (varSum[j] - varSum[j - 1]) / stepsCurrent - \
                            ((xiSum[j] - xiSum[j - 1]) / stepsCurrent) ** 2
                xiSep[j] = (xiSum[j] - xiSum[j - 1]) / stepsCurrent
        # print(varSep)
        # print(xiSep)
        xiSep = np.add(xiSep, -xi_list[i])
        xiSep_mean = np.mean(xiSep)
        varSep_mean = np.mean(varSep)
        xiSep_max, xiSep_min = np.max(xiSep[:timeCount_eff]), np.min(xiSep[:timeCount_eff])
        varSep_max, varSep_min = np.max(varSep[:timeCount_eff]), np.min(varSep[:timeCount_eff])
        xiSep_range = xiSep_max - xiSep_min
        varSep_range = varSep_max - varSep_min

        # mean evolution
        me.plot(timeCount, xiMean, c=color[0])
        me.yaxis.set_major_formatter(myFormatter())
        me_var = me.twinx()
        me_var.plot(timeCount, varMean, c=color[1])
        me_var.yaxis.set_major_formatter(formatter)
        me.yaxis.set_major_formatter(myFormatter())
        # me.axhline(y=xiSep_mean, c='blue', ls='--', zorder=10)
        # me_var.axhline(y=varSep_mean, c='red', ls='--', zorder=10)
        me.axhline(y=0.0, c='black', ls='--', lw=1)  # xi_ref
        # axis label
        me.set_xlabel('Time (ps)')
        me.tick_params('y', colors=color[0])
        # me.set_ylabel('$\\xi - \\xi_{{\mathrm{{ref}}}}, \\xi_{{\mathrm{{ref}}}} = {0:.3f}$'.format(xi_list[i]))
        me.set_ylabel('$\\xi - \\xi_{{\mathrm{{ref}}}}$', color=color[0])
        me_var.tick_params('y', colors=color[1])
        me_var.set_ylabel('Variance', color=color[1])
        # axis lim
        me.set_ylim(xiSep_min - xiSep_range * 0.1, xiSep_max + xiSep_range * 0.1)
        me_var.set_ylim(varSep_min - varSep_range * 0.1, varSep_max + varSep_range * 0.1)

        # mean scatter
        ms.plot(xiMean, varMean, c=color[0], marker='o', ls='', markersize=1)
        ms.yaxis.set_major_formatter(myFormatter())
        ms.xaxis.set_major_formatter(myFormatter())
        ms.axvline(x=0.0, c='black', ls='--', lw=1)  # xi_ref
        ms.set_ylabel('Variance')
        ms.set_xlabel('$\\xi - \\xi_{{\mathrm{{ref}}}}$')
        # axis lim
        ms.set_ylim(varSep_min - varSep_range * 0.1, varSep_max + varSep_range * 0.1)
        ms.set_xlim(xiSep_min - xiSep_range * 0.1, xiSep_max + xiSep_range * 0.1)

        # individual evolution
        ee.plot(timeCount, xiSep, c=color[0], lw=0.5)
        ee.yaxis.set_major_formatter(myFormatter())
        ee_var = ee.twinx()
        ee_var.plot(timeCount, varSep, c=color[1], lw=0.5)
        ee_var.yaxis.set_major_formatter(myFormatter())
        # ee.axhline(y=xiSep_mean, c='blue', ls='--', zorder=10)
        # ee_var.axhline(y=varSep_mean, c='red', ls='--', zorder=10)
        ee.axhline(y=0.0, c='black', ls='--', lw=1)  # xi_ref
        ee.set_xlabel('Time (ps)')
        ee.tick_params('y', colors=color[0])
        ee.set_ylabel('$\\xi - \\xi_{\mathrm{ref}}$', color=color[0])
        ee_var.tick_params('y', colors=color[1])
        ee_var.set_ylabel('Variance', color=color[1])
        # axis lim
        ee.set_ylim(xiSep_min - xiSep_range * 0.1, xiSep_max + xiSep_range * 0.1)
        ee_var.set_ylim(varSep_min - varSep_range * 0.1, varSep_max + varSep_range * 0.1)

        # individual scatter
        es.plot(xiSep, varSep, c=color[0], marker='o', ls='', markersize=1)
        es.xaxis.set_major_formatter(myFormatter())
        es.yaxis.set_major_formatter(myFormatter())
        # es.axhline(y=varSep_mean, c='red', ls='--', zorder=10)
        # es.axvline(x=xiSep_mean, c='blue', ls='--', zorder=10)
        es.axvline(x=0.0, c='black', ls='--', lw=1)  # xi_ref
        es.set_ylabel('Variance')
        es.set_xlabel('$\\xi - \\xi_{{\mathrm{{ref}}}}$')
        # axis lim
        es.set_ylim(varSep_min - varSep_range * 0.1, varSep_max + varSep_range * 0.1)
        es.set_xlim(xiSep_min - xiSep_range * 0.1, xiSep_max + xiSep_range * 0.1)

        # plt.show()
        fig.suptitle(' Variance in window $\\xi_{{\mathrm{{ref}}}}={:.3f}$'.format(xi_list[i]),
                     x=0., y=0.95, horizontalalignment='left', verticalalignment='bottom')
        # fig.suptitle(mylabel,
        #              x=1, y=0.95, horizontalalignment='right', verticalalignment='bottom')
        plot_save('Variances\\{:0>3d}_{:.3f}'.format(i, xi_list[i]))

    return


def plot_deviation():
    title = 'deviation'
    plot_parameters(title)

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    formatter2 = ticker.ScalarFormatter(useMathText=True)
    formatter2.set_scientific(True)
    formatter2.set_powerlimits((-1, 1))

    # NtrajEff
    xi_dev = np.zeros(len(xi_list))
    var = np.zeros(len(xi_list))

    for i in range(len(xi_list)):
        xi_dev[i] = umbInfo[3, NtrajEff - 1, i] - xi_list[i]
        tmp = umbInfo[1, NtrajEff - 1, i] / umbInfo[2, NtrajEff - 1, i] - \
              (umbInfo[0, NtrajEff - 1, i] / umbInfo[2, NtrajEff - 1, i]) ** 2
        var[i] = tmp * kforce_list[i] * 627.509474063056  # to kcal/mol
        # * temp *, no temperature here.

    plt.plot(xi_list, xi_dev, c=color[0])
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.xlabel('$\\xi$')
    plt.ylabel('$\\xi_i - \\xi^{\\mathrm{ref}}_i$')

    plot_var = plt.twinx()
    plot_var.plot(xi_list[0], var[0], c=color[0], label='$\\xi_i - \\xi^{\\mathrm{ref}}_i$')
    plot_var.plot(xi_list, var, c=color[1], label='$\\sigma_i k_i$')
    plot_var.yaxis.set_major_formatter(formatter2)
    plot_var.set_ylabel('$\\sigma_i k_i $ (kcal/mol)')
    # plot_var.set_minorticks_on()

    plot_var.legend(loc='best')
    plot_save('xi_dev')

def plot_pmf(path):
    title = 'PMF'
    plot_parameters(title)

    try:
        f = open(path + '/potential_of_mean_force.dat', 'r')
    except FileNotFoundError:
        print('[ERROR] {} file not found! '.format(title))
    else:
        fLines = f.readlines()
        f.close()

        xi = []
        pmf = []
        for i in fLines[12:-1]:
            xi.append(np.float(i.split()[0]))
            pmf.append(np.float(i.split()[1]))

        # Let W(xi=0) = 0!
        xiAbs = np.abs(xi)
        xiZeroIndex = list(xiAbs).index(min(np.abs(xi)))
        pmf = [(x - pmf[xiZeroIndex]) / 27.211386245988 * 627.509474063056 for x in
               pmf]  # shift and convert to kcal/mol

        # write PMF in kcal/mol
        with open(os.path.join(figPath, 'PMF.txt'), 'w') as pmfFile:
            for i in range(len(xi)):
                pmfFile.write('{:.6f}\t{:.10f}\n'.format(xi[i], pmf[i]))

        plt.xlabel(r'Reaction Coordinate')
        plt.ylabel(r'$W(\xi)$ (kcal/mol)')

        # Choose the fit range of this plot
        plt.xlim(xi[0], xi[-1])
        yRange = max(pmf) - min(pmf)
        plt.ylim(min(pmf) - yRange * 0.1,
                 max(pmf) + yRange * 0.1)  # adds 0.1*yRange to the top and bottom

        plt.plot(xi, pmf, c=color[0], label=mylabel)

        # # plot a zoomed subfigure
        # xiMaxIndex = pmf.index(max(pmf)) # the position of maximum
        # extend = 80 # Extra points to plot
        # ximax = xi[xiMaxIndex - extend : xiMaxIndex + extend]
        # pmfmax = pmf[xiMaxIndex - extend : xiMaxIndex + extend]
        #
        # # # Find maximum of xi
        # f = itp(ximax, pmfmax, k=4)
        # cr_pts = f.derivative().roots()
        # cr_pts = np.append(cr_pts, (ximax[0], ximax[-1]))  # also check the endpoints of the interval
        # cr_vals = f(cr_pts)
        # #! min_index = np.argmin(cr_vals)
        # max_index = np.argmax(cr_vals)
        # pmfMax, xiMax = cr_vals[max_index], cr_pts[max_index]
        #
        # subfig = plt.axes([.3, .5, .5, .4])
        # subfig.plot(ximax, pmfmax, c=color[0])
        #
        # subfig.axvline(x=xiMax, c=color[0], lw=0.5, linestyle='--')
        # subfig.axhline(y=pmfMax, c=color[0], lw=0.5, linestyle='--')
        #
        # plt.setp(subfig, xlim=[min(ximax), max(ximax)])

        plt.legend(loc='upper left')
        plot_save(title)


def plot_rexFactor(path):
    title = 'Transmission_Coefficient'
    plot_parameters(title)

    # Find the file first!
    fileList = os.listdir(path)
    rexFileName = ''
    for file in fileList:
        if file[:18] == 'recrossing_factor_':
            rexFileName = file
    if rexFileName=="":
        pass
    else:
        try:
            f = open(path + '/' + rexFileName, 'r')
        except FileNotFoundError:
            print('[ERROR] {} file not found! '.format(title))
        except PermissionError:
            print('[ERROR] {} file not found! '.format(title))
        else:
            fLines = f.readlines()
            f.close()
            time = []
            kappa = []
            for i in fLines[17:-1]:
                ele = i.split()
                time.append(np.float(ele[0]))
                kappa.append(np.float(ele[-1]))

            plt.xlabel('$t$ (fs)')
            plt.ylabel('$\kappa(t)$')

            # plt.xscale('log')

            plt.xlim(time[0], time[-1])
            # # endRF = np.mean(kappa[-5:])
            # plt.axhline(y=kappa[-1], c=color[0], lw=0.5, linestyle='--')

            plt.plot(time, kappa, c=color[0], label=mylabel)

            plt.legend(loc="best")
            plot_save(title)

            with open(os.path.join(figPath, 'recrossing.txt'), 'w') as rexFile:
                for i in range(len(time)):
                    rexFile.write('{:.3f}\t{:.6f}\n'.format(time[i], kappa[i]))


def plot_overlap_density(path):
    if NtrajEff == 1:
        return

    title = 'Overlap_Density'
    plot_parameters(title)

    resolution = 2000
    extend = 0.03  # 3E-2

    xiMin = np.min(xi_list)
    xiMax = np.max(xi_list)

    sizeV = NtrajEff  # np.shape(umbInfo)[1]
    sizeH = np.shape(umbInfo)[2]

    x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)

    # Read time unit
    tempFile = open(path + "\\mbrella_sampling_{0:.8f}.dat".format(xi_list[0]), 'r')
    lines = tempFile.readlines()
    timeSep = np.float(lines[9].split()[4]) / 1000.0  # to ns

    z = np.zeros((sizeV * resolution))
    y = np.linspace(xiMin - extend, xiMax + extend, resolution)

    # Gaussian summation
    for j in range(sizeV):  # xi
        y_sum = np.zeros((resolution))
        for i in range(sizeH):  # var
            y_new = my_gaussian(x_new, umbInfo[3, j, i], umbInfo[4, j, i])
            y_sum += y_new
            z[j * resolution:(j + 1) * (resolution)] = y_sum


    # plt.show()
    z = z.reshape(sizeV, resolution).transpose()
    z /= np.max(z)
    x = np.linspace(0, timeSep * sizeV, sizeV)
    plt.pcolormesh(x, y, z, cmap='Greens', vmax=1.0)  # pcolormesh # contourf # Greys_r

    plt.xlabel('Time (ns)')
    plt.ylabel('Reaction Coordinate')

    plt.colorbar()
    # plt.title('The Evolution of Normalized Population')
    plot_save(title)

    # 3D UI
    plot_parameters('UI (3D)')
    X, Y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(5, 3.75))  # 1.25 * (4,3)
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, z, cmap='Greens', linewidth=0.2, edgecolors='black')
    ax.view_init(elev=20, azim=30)

    # ax.set_zticks([])
    ax.set_xlabel(r'Time (ns)')
    ax.set_ylabel(r'Reaction Coordinate')
    ax.set_zlabel(r'Normalized Population')

    plot_save('Overlap_Density_3D')

    clearFolder('UI')
    for cycle in range(NtrajEff):
        if (cycle + 1) % np.ceil(NtrajEff / 10) == 0 or cycle == 0 or cycle == NtrajEff - 1:
            resolution = 2000
            extend = 0.03  # 3E-2

            xiMin = np.min(xi_list)
            xiMax = np.max(xi_list)
            length = len(xi_list)

            x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)
            y_sum = np.zeros((resolution))  # Total density line

            maxPop = 0  # maximum of the summation of all population

            timeCurrent = umbInfo[2, cycle, 0] * delta  # * 1E-3 # to ps
            plot_parameters('UI at time {:.4f} ps'.format(timeCurrent))
            plt.figure(figsize=(9, 3))

            for i in range(length):
                # Gaussian smearing
                xav, xav2 = umbInfo[3:, cycle, i]
                y_new = my_gaussian(x_new, xav, xav2)

                if xav2 - 1E-8 < 0:
                    print(cycle, i, '<0', xav, xav2)

                # Find biggest population
                if max(y_new) > maxPop:
                    maxPop = max(y_new)

                y_sum += y_new  # sum all population
                if xav2 > 5.0E-5:
                    # print("[WARNING] May be too various in xi = {}! ".format(xi_list[i]))
                    plt.plot(x_new, y_new, lw=1, c=color[1], alpha=0.8)
                else:
                    plt.plot(x_new, y_new, lw=0.5, c=color[0], alpha=.3)

            # Plot summation and difference
            plt.plot(x_new, y_sum, lw=1, c=color[0], label='{:.0f} ps'.format(timeCurrent))
            # mylabel)  # label='Summation of all populations')  # SHU Blue

            # plt.xlabel('Reaction Coordinate / Å')
            plt.xlabel('Reaction Coordinate')
            plt.ylabel('Population')

            plt.xlim(xiMin - extend, xiMax + extend)
            # plt.ylim(0, maxPop*1.1)
            plt.ylim(0, max(y_sum) * 1.2)

            plt.yticks([])  # No ticks and labels in y axis

            plt.legend(loc='upper left')
            plot_save('UI\\{:.0f}'.format(timeCurrent))



def plotKForce():
    plot_parameters('force constant')

    markerline, stemlines, baseline = \
        plt.stem(xi_list, kforce_list, use_line_collection=True,
                 basefmt=' ', markerfmt=' ', linefmt=color[0])
    plt.setp(stemlines, 'linewidth', 0.5)
    plt.scatter(xi_list, kforce_list, c=color[0], s=1, label=mylabel)

    # upper bound of force constants
    kMax = np.zeros(len(xi_list))
    for i in range(len(xi_list) - 1):
        kMax[i] = 9 * 3.16681046247368 * 1E-6 / (xi_list[i + 1] - xi_list[i]) ** 2
        # print(kMax[i])
    kMax[-1] = kMax[-2]
    plt.plot(xi_list, kMax, '--', lw=0.75, c=color[1], alpha=0.5)
    # emphasis the large force constants
    for i, kmax in enumerate(kMax):
        if kforce_list[i] > kmax:
            print('       Large kforce {:.3f} in xi={:.3f}'.format(kforce_list[i], xi_list[i]))
            plt.scatter(xi_list[i], kforce_list[i], c='red', s=4, zorder=10)
            # markerline, stemlines, baseline = \
            #     plt.stem(xi_list[i], kforce_list[i], use_line_collection=True,
            #              basefmt=' ', markerfmt=' ', linefmt='red')
            # plt.setp(stemlines, 'linewidth', 0.5)

    plt.ylabel('Force Constant (Hartree)')  # $T$ K$^{-1}$
    plt.xlabel('Reaction Coordinate')

    # plt.xlim(min(xi_list), max(xi_list))
    plt.ylim(0, max(kforce_list) * 1.1)

    # plt.legend(loc='lower right')
    plt.legend(loc='upper left')
    plot_save('kforce')


def plot_PMF_evolution():
    if NtrajEff == 1:
        return

    print('[INFO] Computing PMF evolution...')
    clearFolder('PMF')

    # Constants
    bins = 500
    beta = 4.35974417e-18 / (1.3806504e-23 * temp)
    totalCycle = NtrajEff  #np.shape(umbInfo)[1]  # the number of trajectories
    Nwindows = np.shape(umbInfo)[2]  # the number of windows

    # middle variables for calculations
    binList = np.linspace(min(xi_list), max(xi_list), bins, True)
    dA = np.zeros(bins)
    p = np.zeros(Nwindows)  # probability
    dA0 = np.zeros(Nwindows)

    # PMF data storage
    PMFdata = np.zeros((bins - 1, totalCycle))
    freeEnergy = np.zeros((3, totalCycle))  # time, xi, free energy

    for cycle in range(totalCycle):
        # save 10 PMF figures
        if True:#(cycle + 1) % np.ceil(totalCycle / 10) == 0 or cycle == 0 or cycle == totalCycle:
            # for cycle in [-1]:
            print('       Computing PMF evolution {} of {}'.format(cycle + 1, totalCycle))
            N = umbInfo[2, cycle, :]
            for n, xi in enumerate(binList):
                for l in range(Nwindows):
                    # av = window.av / window.count
                    # av2 = window.av2 / window.count
                    xi_mean = umbInfo[3, cycle, l]  # computed
                    xi_var = umbInfo[4, cycle, l]  # computed
                    xi_window = xi_list[l]  # fixed
                    kforce = kforce_list[l] * temp  # fixed
                    p[l] = 1.0 / np.sqrt(2 * np.pi * xi_var) * np.exp(-0.5 * (xi - xi_mean) ** 2 / xi_var)
                    dA0[l] = (1.0 / beta) * (xi - xi_mean) / xi_var - kforce * (xi - xi_window)
                    # plt.plot(p, label=l)
                dA[n] = np.sum(N * p * dA0) / np.sum(N * p)
                # plt.legend()
                # plt.show()

            # Now integrate numerically to get the potential of mean force
            potentialOfMeanForce = np.zeros((2, bins - 1))
            A = 0.0
            for n in range(bins - 1):
                dx = binList[n + 1] - binList[n]
                potentialOfMeanForce[0, n] = 0.5 * (binList[n] + binList[n + 1])
                A += 0.5 * dx * (dA[n] + dA[n + 1])
                potentialOfMeanForce[1, n] = A
            potentialOfMeanForce[1, :] -= np.min(potentialOfMeanForce[1, :])

            PMFcurrent = potentialOfMeanForce[1, :] * 627.503  # to kcal/mol

            # Let W(xi=0) = 0!
            xiAbs = np.abs(binList)
            xiZeroIndex = list(xiAbs).index(min(np.abs(binList)))
            PMFcurrent = [x - PMFcurrent[xiZeroIndex] for x in PMFcurrent]
            PMFdata[:, cycle] = PMFcurrent

            timeCurrent = umbInfo[2, cycle, 0] * delta  # * 1E-3


            plot_parameters('PMF at time {:.0f} ps'.format(timeCurrent))
            plt.plot(binList[:-1], PMFcurrent, c=color[0], label='{:.0f} ps'.format(timeCurrent))
            plt.xlabel(r'Reaction Coordinate')
            plt.ylabel(r'$W(\xi)$ (kcal/mol)')
            plt.legend(loc='upper left')
            plot_save('PMF\\{:.0f}'.format(timeCurrent))

        # calculate free energy
        pmfMaxValue = np.max(PMFcurrent)
        pmfMaxIndex = PMFcurrent.index(pmfMaxValue)
        freeEnergy[:, cycle] = timeCurrent, binList[pmfMaxIndex], pmfMaxValue


    # # write PMF datas
    # f = open('PMF_data.txt', 'w')
    # f.write('\t'.join(map(str, list(umbInfo[2,:,0]))) + '\n') # map(str, value_list)
    # for i in range(bins - 1):
    #     f.write('\t'.join(map(str, PMFdata[i,:])) + '\n')

    # write PMF datas
    f = open(os.path.join(figPath, 'PMF_data.txt'), 'w')
    for j in range(totalCycle):
        for i in range(bins - 1):
            f.write(
                '{:.4f}\t{:.4f}\t{:.4f}\n'.format(umbInfo[2, j, 0] * delta, binList[i],  # * 1E-3
                                                  PMFdata[
                                                      i, j]))  # time to ns # ps # 2020-05-02 15:44:43 Wenbin, FAN @ SHU
    f.close()

    # Plot PMF evolution
    plot_parameters('PMF evolution')

    a = pd.read_csv(os.path.join(figPath, 'PMF_data.txt'), sep='\t', header=None)

    traj = list(a.iloc[:, 0])
    xibins = list(a.iloc[:, 1])
    pmfvalue = list(a.iloc[:, 2])

    pmfMin = np.min(pmfvalue)
    pmfMax = np.max(pmfvalue)
    level = np.arange(int(pmfMin) - 2, int(pmfMax) + 2, 1)

    plt.tricontourf(traj, xibins, pmfvalue, levels=level, cmap='Blues')
    plt.colorbar()
    plt.tricontour(traj, xibins, pmfvalue, linestyles='-', levels=level, colors='Black', linewidths=0.2)

    plt.xlabel(r'Time (ps)')
    plt.ylabel(r'Reaction Coordinate')
    plot_save('PMF_evolution')

    # 3D plot
    # print(len(traj), (totalCycle, bins - 1))
    X = np.reshape(traj, (totalCycle, bins - 1))
    Y = np.reshape(xibins, (totalCycle, bins - 1))
    Z = np.reshape(pmfvalue, (totalCycle, bins - 1))

    fig = plt.figure(figsize=(5, 3.75))  # 1.25 * (4,3)
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, cmap='Blues', linewidth=0.2, edgecolors='black')
    ax.view_init(elev=20, azim=30)

    ax.set_xlabel(r'Time (ps)')
    ax.set_ylabel(r'Reaction Coordinate')
    ax.set_zlabel(r'Free Energy (kcal/mol)')

    plot_save('PMF_evolution_3D')

    # Plot free energy
    plot_parameters('free energy')

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    fig.set_figheight(3)
    fig.set_figwidth(4)
    ax1.set_xlim(0, np.max(freeEnergy[0, :]))

    ax1.plot(freeEnergy[0, :], freeEnergy[1, :], c=color[0])  # , label='Reaction Coordinate')
    ax2.plot(freeEnergy[0, :], freeEnergy[2, :], c=color[1], label='Free Energy')
    ax2.plot(freeEnergy[0, 0], freeEnergy[2, 0], c=color[0], label='Reaction Coordinate')  # fake figure for legend

    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Reaction Coordinate')  # , color=color[0])
    ax2.set_ylabel('Free Energy (kcal/mol)')  # , color=color[1])

    ax2.legend(loc='best')

    plot_save('PMF_free_energy')


def plot_xi():
    plot_parameters('xi evolution')

    length = len(xi_list)
    xiref_evolution = np.zeros((Ntraj, length))

    timeMax = 0.0
    timeMin = 1E5
    for i in range(length):
        tscolor = (((Tcolor2[0] - Tcolor1[0]) * i / length + 0) / 255.0,
                   ((Tcolor2[1] - Tcolor1[1]) * i / length + 68) / 255.0,
                   ((Tcolor2[2] - Tcolor1[2]) * i / length + 124) / 255.0)

        timeEvolution = umbInfo[2, :, i]
        timeEvolution = [x * delta for x in
                         timeEvolution]  # 0.1 fs to 1 ns #  * 1E-3 # ps # 2020-05-02 15:46:21 Wenbin, FAN @ SHU
        xiEvolution = umbInfo[3, :, i]
        xiref_evolution[:, i] = np.add(np.divide(umbInfo[0, :, i], umbInfo[2, :, i]), -xi_list[i])

        # light color for normal xi
        alpha = 0.0
        if max(xiEvolution) - min(xiEvolution) > (xi_list[1] - xi_list[0]) / 5.0:
            alpha = 1.0
        else:
            alpha = 0.3

        plt.plot(timeEvolution, umbInfo[3, :, i], c=tscolor, lw=0.5, alpha=alpha)

        if max(timeEvolution) > timeMax:
            timeMax = max(timeEvolution)
        if min(timeEvolution) < timeMin:
            timeMin = min(timeEvolution)

    plt.xlim(0, timeMax)  # timeMin, timeMax
    plt.xlabel('Time (ps)')
    plt.ylabel('Reaction Coordinates')

    plot_save('xi_evolution')

    plot_parameters('xi-ref_evolution')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    for i in range(length):
        tscolor = (int((255. * i / length)) / 255.0,
                   0.,
                   (int(-255. * i / length + 255)) / 255.0)
        plt.plot(timeEvolution, xiref_evolution[:, i], c=tscolor)

    plt.xlim(0, timeMax)
    plt.xlabel('Time (ps)')
    plt.ylabel('$\\xi_i - \\xi_i^{\\mathrm{ref}}$')
    plt.gca().yaxis.set_major_formatter(formatter)
    plot_save('xi-ref_evolution')

    return


def getBasicInfo(path):
    global args
    T=args.T
    N=args.N
    
    """# get the name of submitting script
    subList = ['run.sh', 'highcpu', 'fat', 'gpu', 'pbs', 'run.txt', 'sub.lsf', 'sub.pbs', 'fat.yy']  # submitting script
    subName = ''
    subPath = ''
    for i, name in enumerate(subList):
        if os.path.exists(os.path.join(path, name)):
            subName = name
            subPath = os.path.join(path, name)
            break

    # read temperature and the number of beads
    try:
        f = open(subPath, 'r')
    except FileNotFoundError:
        f = open("input.py", 'r')

    print('[INFO] Submitting arguments: ')
    cmdLine = ''
    for line in f:
        if line[:6] == 'python':
            cmdLine = line.split()
    del cmdLine[:2]  # delete `python` and `...rpmdrate.py`

    # get the number of cores used
    cores = 0
    for i, cmd in enumerate(cmdLine):
        if cmd == '-p':
            cores = cmdLine[i + 1]
            print('       Cores:            {}'.format(int(cores)))
            del cmdLine[i:i + 2]
            break
    if cores == 0:
        print('       Cores:            single')

    # delete redirect output
    for i, cmd in enumerate(cmdLine):
        if cmd == '>' or cmd == '>>' or cmd == '#>' or cmd == '|':
            del cmdLine[i:]

    # get input file, temperature and the number of bsubmitting scripteads
    assert len(cmdLine) == 3, 'Your submitting script may be wrong! '"""
    global inputFile, temp, Nbeads  # There will be probably more pythonic way. Tell me plz if you know!
    inputFile = args.I
    temp = T
    Nbeads = N

    print('       Temperature:      {} K'.format(temp))
    print('       Number of beads:  {}'.format(Nbeads))
    print('       Input file:       {}\n'.format(inputFile))

    return inputFile


def getUmbrellaInfo(path):
    print('[INFO] Getting umbrella data...')
    Nwindows = len(xi_list)
    print('       number of windows: {}'.format(Nwindows))

    # Count the total lines of xi and xvar
    NtrajList = np.zeros(Nwindows)
    for i in range(Nwindows):
        with open(path + "/umbrella_sampling_{0:.8f}.dat".format(xi_list[i]), 'r') as tempFile:
            for j, l in enumerate(tempFile):
                pass
        NtrajList[i] = j + 1 - 15  # 15 info lines # +1 means the number of lines

    global Ntraj, NtrajEff
    Ntraj = int(np.max(NtrajList))
    NtrajEff = int(np.min(NtrajList))  # Effective lines
    print('       Maximum of trajectories: {}'.format(Ntraj))
    print('       Minimum of trajectories: {}\n'.format(NtrajEff))

    global umbInfo
    umbInfo = np.zeros((5, Ntraj, Nwindows))  # `5` means five columns in the umbrella info files.

    # Read time unit
    tempFile = open(path + "/umbrella_sampling_{0:.8f}.dat".format(xi_list[0]), 'r')
    lines = tempFile.readlines()
    timeSep = np.float(lines[9].split()[4])  # / 1000.0  # to ns # ps # 2020-05-02 15:46:42 Wenbin, FAN @ SHU
    tempFile.close()

    # Read in all data
    for i in range(Nwindows):
        fname = path + "/umbrella_sampling_{0:.8f}.dat".format(xi_list[i])
        f = open(fname, 'r').readlines()

        for j in range(Ntraj):
            try:
                line = f[15 + j].split()
                if len(line) != 5:
                    print(f'{xi_list[i]} has wrong line! ')
                    raise Exception
                umbInfo[:, j, i] = line
            except IndexError:
                umbInfo[:, j, i] = None

    return umbInfo


def getInput(folder):
    inputPath =inputFile

    # skip `import PES`
    f = open(inputPath, 'r')
    start = 0
    inputContent = ''
    for i, line in enumerate(f.readlines()):
        if line[:5] == 'label' or line[:5] == 'react':
            start = 1
        if start == 1:
            inputContent += line.replace('numpy', 'np')

    # execute the input.py and get force constant
    T = temp
    try:
        exec(inputContent)
    except (NameError, TypeError, SyntaxError):
        print('[ERROR] The input file {0!r} was invalid:'.format(inputPath))
        raise

    global path
    path = os.path.join(folder, str(temp), str(Nbeads))

    global mylabel
    if Nbeads == 1:
        mylabel = '{} K, {} bead'.format(temp, Nbeads)
    else:
        mylabel = '{} K, {} beads'.format(temp, Nbeads)


def getRate(path):
    print('[INFO] rate coefficients: ')
    fileList = os.listdir(path)
    rateFile = ''
    for file in fileList:
        if file.split('_')[0] == 'rate':
            rateFile = os.path.join(path, file)
            break

    if len(rateFile) == 0:
        print('[INFO] No rate file. ')
        return

    f = open(rateFile, 'r')
    g = open(os.path.join(figPath, 'my_rate.txt'), 'w')
    fl = f.readlines()

    rateTemp = np.float(fl[4].split()[-2])
    rateProb = np.float(fl[10].split()[-1])
    rateMaxxi = np.float(fl[12].split()[-1])
    rateQTST = np.float(fl[14].split()[-2])
    rateRex = np.float(fl[17].split()[-1])
    rateRPMD = np.float(fl[19].split()[-2])

    rateFreeEnergy = - np.log(rateProb) * rateTemp * 1.3806504E-23 / 4.3597447222071E-18  # Hartree
    # # kB = 1.3806504E-23 (J/K), 1 Hartree = 4.3597447222071e-18 (J)
    rateRPMDfT = rateRPMD * 2 / (2 + 2 * np.exp(- 205 / rateTemp))

    print('Temperature (K):   \t{:d}'.format(int(rateTemp)))
    print('xi^ddagger:        \t{:.3f}'.format(rateMaxxi))
    print('delta G (kcal/mol):\t{:.2f}'.format(rateFreeEnergy * 627.509474063056))
    print('k_QTST:            \t{:.2e}'.format(rateQTST))
    print('kappa:             \t{:.3f}'.format(rateRex))
    print('k_RPMD:            \t{:.2e}'.format(rateRPMD))
    print('k_RPMD * f(T):     \t{:.2e}'.format(rateRPMDfT))

    g.write('Temperature (K):   \t{:d}  \n'.format(int(rateTemp)))
    g.write('xi^ddagger:        \t{:.3f}\n'.format(rateMaxxi))
    g.write('delta G (kcal/mol):\t{:.2f}\n'.format(rateFreeEnergy * 627.509474063056))
    g.write('k_QTST:            \t{:.2e}\n'.format(rateQTST))
    g.write('kappa:             \t{:.3f}\n'.format(rateRex))
    g.write('k_RPMD:            \t{:.2e}\n'.format(rateRPMD))
    g.write('k_RPMD * f(T):     \t{:.2e}\n'.format(rateRPMDfT))

    f.close()
    g.close()

    return

# Defination in RPMDrate:
def reactants(atoms, reactant1Atoms, reactant2Atoms, Rinf):
    pass


def transitionState(geometry, formingBonds, breakingBonds):
    pass


def equivalentTransitionState(formingBonds, breakingBonds):
    pass


def thermostat(type, **kwargs):
    pass


def generateUmbrellaConfigurations(dt, evolutionTime, xi_list, kforce):
    global delta
    if dt[1] == 'ps':
        delta = dt[0]
    elif dt[1] == 'fs':
        delta = dt[0] * 1E-3
    else:
        print('[ERROR] Time unit {} not support and will be regarded as `ps`. '.format(dt[1]))
        delta = dt[0]
    assert np.float(delta) < 1

def conductUmbrellaSampling(dt, windows, saveTrajectories=False):
    global xi_list, kforce_list
    # print(windows)
    xi_list = np.zeros(len(windows))
    kforce_list = np.zeros(len(windows))
    for i in range(len(windows)):
        xi_list[i] = '{0:.8f}'.format(windows[i][0])
        kforce_list[i] = windows[i][1] / temp

    kf_path = os.path.join(os.path.abspath(os.path.dirname(figPath)), 'kforce.txt')
    if os.path.exists(kf_path):
        kf_list_read = open(kf_path, 'r')
        kflines = kf_list_read.readlines()
        for i, line in enumerate(kflines):
            kforce_list[i] = np.float(line.split()[1])
            # print('kforce: {}'.format(kforce_list[i]))
    return xi_list, kforce_list

def computePotentialOfMeanForce(windows=None, xi_min=None, xi_max=None, bins=5000):
    pass


def computeRecrossingFactor(dt, equilibrationTime, childTrajectories, childSamplingTime, childrenPerSampling,
                            childEvolutionTime, xi_current=None, saveParentTrajectory=False,
                            saveChildTrajectories=False):
    pass


def computeRateCoefficient():
    pass


def Window(xi, kforce, trajectories, equilibrationTime, evolutionTime):
    return [np.float(xi), np.float(kforce)]


def main(inputFolder=None):
    # if inputFolder == None:
    #     inputFolder = input_path()
    inputFolder=args.R
    if inputFolder[-1]!="/":
        inputFolder+="/"
    global figPath
    figPath = os.path.join(inputFolder, str(args.T)+"_"+str(args.N)+'fig')
    if not os.path.exists(figPath):
        os.mkdir(figPath)

    # get info
    getBasicInfo(inputFolder)
    getInput(inputFolder)
    getUmbrellaInfo(path)
    getRate(path)

    # # plot
    plotKForce()
    plot_overlap()
    plot_variance()
    plot_variance_diff()
    plot_pmf(path)
    plot_rexFactor(path)
    plot_xi()
    plot_deviation()

    # plot_PMF_evolution()
    # plot_var_evolution()
    # plot_overlap_density(path)

    myEnding()

main()

# root = r'C:\Users\Mike\Desktop\fin-OD-300_2'
# for dir in os.listdir(root):
#     print(dir)
#     main(os.path.join(root, dir))
