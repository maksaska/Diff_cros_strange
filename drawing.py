from model import Model
import os.path
from math import *
import numpy as np
import matplotlib.pyplot as plt

#
#
# PLOTTING FUNCTIONS
#
#


def plot_cs_W(obj: Model, Q2: float, cos_th: float, E_beam: float, phi: float):
    obj.Q2 = Q2
    obj.cos_th = cos_th
    obj.E_beam = E_beam
    obj.phi = phi

    if obj.name == "K+L":
        Wmin = 1.61
    else:
        Wmin = 1.69

    x = np.arange(Wmin, 2.65, 0.01)
    y = []
    dy = []

    for W in x:
        obj.W = W
        res = obj.Point_diff()
        y.append(res[0])
        dy.append(res[1])

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'$Q^{2}$ = {Q2} GeV$^{2}$\n $\cos$'+r'$\theta$ = ' + f'{cos_th}\nE = {E_beam} GeV\n $\phi$ = {phi} degree')

    plt.xlabel('W, GeV')
    plt.ylabel('$d\sigma/d\Omega$, [nb/sr]')
    plt.title('W - dependence')
    plt.legend(prop={'size': 20})
    plt.plot([Wmin, 2.65], [0, 0], color='black')

    # plt.show()

    plt.savefig(
        f"./pictures/cs_W_{obj.name}_{Q2}_{cos_th}_{E_beam}_{phi}.pdf", bbox_inches='tight')
    plt.clf()


def plot_cs_Q2(obj: Model, W: float, cos_th: float, E_beam: float, phi: float):
    obj.W = W
    obj.cos_th = cos_th
    obj.E_beam = E_beam
    obj.phi = phi

    x = np.arange(0.0, 5.0, 0.05)
    y = []
    dy = []

    for Q2 in x:
        obj.Q2 = Q2
        res = obj.Point_diff()
        y.append(res[0])
        dy.append(res[1])

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'W = {W} GeV\n $\cos$'+r'$\theta$ = ' + f'{cos_th}\nE = {E_beam} GeV\n $\phi$ = {phi} degree')

    plt.xlabel('$Q^{2}$, GeV$^{2}$')
    plt.ylabel('$d\sigma/d\Omega$, [nb/sr]')
    plt.title('$Q^{2}$ - dependence')
    plt.legend(prop={'size': 20})
    plt.plot([0, 5.0], [0, 0], color='black')

    # plt.show()

    plt.savefig(
        f"./pictures/cs_Q2_{obj.name}_{W}_{cos_th}_{E_beam}_{phi}.pdf", bbox_inches='tight')
    plt.clf()


def plot_cs_cos(obj: Model, W: float, Q2: float, E_beam: float, phi: float):
    obj.W = W
    obj.Q2 = Q2
    obj.E_beam = E_beam
    obj.phi = phi

    x = np.arange(-1, 1, 0.01)
    y = []
    dy = []

    for cos_th in x:
        obj.cos_th = cos_th
        res = obj.Point_diff()
        y.append(res[0])
        dy.append(res[1])

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'W = {W} GeV\n$Q^{2}$ = {Q2} GeV$^{2}$ \nE = {E_beam} GeV\n $\phi$ = {phi} degree')

    plt.xlabel(r'$\cos \theta$')
    plt.ylabel('$d\sigma/d\Omega$, [nb/sr]')
    plt.title(r'$\cos \theta$ - dependence')
    plt.legend(prop={'size': 20})
    plt.plot([-1, 1], [0, 0], color='black')

    # plt.show()

    plt.savefig(
        f"./pictures/cs_cos_{obj.name}_{W}_{Q2}_{E_beam}_{phi}.pdf", bbox_inches='tight')
    plt.clf()


def plot_cs_phi(obj: Model, W: float, Q2: float, cos_th: float, E_beam: float):
    obj.W = W
    obj.Q2 = Q2
    obj.E_beam = E_beam
    obj.cos_th = cos_th

    x = obj.phi
    y = []
    dy = []

    res = obj.Point_diff()

    for i in res:
        y.append(i["f"])
        dy.append(i["df"])

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'W = {W} GeV\n$Q^{2}$ = {Q2} GeV$^{2}$ \nE = {E_beam} GeV\n $\cos$' + r'$\theta$ = ' + f'{cos_th}')

    plt.xlabel('$\phi$, degree')
    plt.ylabel('$d\sigma/d\Omega$, [nb/sr]')
    plt.title('$\phi$ - dependence')
    plt.legend(prop={'size': 20})
    plt.plot([-180, 180], [0, 0], color='black')

    # plt.show()

    plt.savefig(
        f"./pictures/cs_phi_{obj.name}_{W}_{Q2}_{E_beam}_{cos_th}.pdf", bbox_inches='tight')
    plt.clf()


def plot_str_W(name: str, obj: Model, Q2: float, cos_th: float, E_beam: float):

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    if name == "St":
        a = 1
    if name == "Sl":
        b = 1
    if name == "Slt":
        c = 1
        j = 6
    if name == "Stt":
        d = 1
        j = 8
    if name == "Su":
        e = 1
        j = 4

    obj.Q2 = Q2
    obj.cos_th = cos_th
    obj.E_beam = E_beam

    if obj.name == "K+L":
        Wmin = 1.61
    else:
        Wmin = 1.69

    x = np.arange(Wmin, 2.65, 0.01)
    y = []
    dy = []

    for W in x:
        obj.W = W
        res = obj.Str_func_all()
        y.append(a*res[0] + b*res[2] + c*res[4] + d*res[6] + e*(res[0] + obj.eps(W, Q2)*res[2]))
        dy.append(a*res[1] + b*res[3] + c*res[5] + d*res[7] + e *
                  sqrt(res[1]**2 + obj.eps(W, Q2)**2*res[3]**2))

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'$Q^{2}$ = {Q2} GeV$^{2}$\n $\cos$'+r'$\theta$ = ' + f'{cos_th}\nE = {E_beam} GeV')

    data_presence = True

    file1 = []
    file2 = []
    file3 = []
    file = []

    if obj.name == "K+L":
        plt.plot(x, y, label=f"$K^+\Lambda$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data_dat/P1.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data_dat/P2.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data_dat/P3.dat', unpack=True)
        if Q2 == 0 and a == 1:
            del file
            file = np.loadtxt('./Data_dat/Diff_L_Photo.dat', unpack=True)
    if obj.name == "K+S0":
        plt.plot(x, y, label=f"$K^+\Sigma^0$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data/P4.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data/P5.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data/P6.dat', unpack=True)
        if Q2 == 0 and a == 1:
            del file
            file = np.loadtxt('./Data/Diff_S_Photo.dat', unpack=True)

    x_data = []
    y_data = []
    dy_data = []

    if a == 1 or b == 1 or len(file1) + len(file2) + len(file3) == 0:
        data_presence = False

    if data_presence:
        if len(file1) > 0:
            for i in range(0, len(file1[1])):
                if file1[1][i] == Q2 and file1[3][i] == cos_th:
                    x_data.append(file1[0][i])
                    y_data.append(file1[j][i])
                    dy_data.append(file1[j+1][i])
        if len(file2) > 0:
            for i in range(0, len(file2[1])):
                if file2[1][i] == Q2 and file2[3][i] == cos_th:
                    x_data.append(file2[0][i])
                    y_data.append(file2[j][i])
                    dy_data.append(file2[j+1][i])
        if len(file3) > 0:
            for i in range(0, len(file3[1])):
                if file3[1][i] == Q2 and file3[3][i] == cos_th:
                    x_data.append(file3[0][i])
                    y_data.append(file3[j][i])
                    dy_data.append(file3[j+1][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS data")

    if a == 1 and Q2 == 0:
        for i in range(0, len(file[1])):
            if abs(file[1][i] - cos_th) < 1e-9:
                x_data.append(file[0][i])
                y_data.append(file[2][i])
                dy_data.append(file[3][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS Photo data")

    plt.xlabel('W, GeV')
    if name == "St":
        plt.ylabel('$S_{t}$, [nb/sr]')
        plt.title('$S_{t}$, [nb/sr]')
    if name == "Sl":
        plt.ylabel('$S_{l}$, [nb/sr]')
        plt.title('$S_{l}$, [nb/sr]')
    if name == "Slt":
        plt.ylabel('$S_{lt}$, [nb/sr]')
        plt.title('$S_{lt}$, [nb/sr]')
    if name == "Stt":
        plt.ylabel('$S_{tt}$, [nb/sr]')
        plt.title('$S_{tt}$, [nb/sr]')
    if name == "Su":
        plt.ylabel('$S_{u}$, [nb/sr]')
        plt.title('$S_{u}$, [nb/sr]')

    plt.legend(prop={'size': 20})

    # plt.show()

    if name == "St":
        plt.savefig(
            f"./pictures/St_W_{obj.name}_{Q2}_{cos_th}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Sl":
        plt.savefig(
            f"./pictures/Sl_W_{obj.name}_{Q2}_{cos_th}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Slt":
        plt.savefig(
            f"./pictures/Slt_W_{obj.name}_{Q2}_{cos_th}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Stt":
        plt.savefig(
            f"./pictures/Stt_W_{obj.name}_{Q2}_{cos_th}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Su":
        plt.savefig(
            f"./pictures/Su_W_{obj.name}_{Q2}_{cos_th}_{E_beam}.pdf", bbox_inches='tight')

    plt.clf()


def plot_str_cos(name: str, obj: Model, W: float, Q2: float, E_beam: float):

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    if name == "St":
        a = 1
    if name == "Sl":
        b = 1
    if name == "Slt":
        c = 1
        j = 6
    if name == "Stt":
        d = 1
        j = 8
    if name == "Su":
        e = 1
        j = 4

    obj.Q2 = Q2
    obj.W = W
    obj.E_beam = E_beam

    x = np.arange(-1, 1, 0.01)
    y = []
    dy = []

    for cos_th in x:
        obj.cos_th = cos_th
        res = obj.Str_func_all()
        y.append(a*res[0] + b*res[2] + c*res[4] + d*res[6] + e*(res[0] + obj.eps(W, Q2)*res[2]))
        dy.append(a*res[1] + b*res[3] + c*res[5] + d*res[7] + e *
                  sqrt(res[1]**2 + obj.eps(W, Q2)**2*res[3]**2))

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'$Q^{2}$ = {Q2} GeV$^{2}$\n W = {W} GeV\nE = {E_beam} GeV')

    data_presence = True

    file1 = []
    file2 = []
    file3 = []
    file = []

    if obj.name == "K+L":
        plt.plot(x, y, label=f"$K^+\Lambda$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data_dat/P1.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data_dat/P2.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data_dat/P3.dat', unpack=True)
        if Q2 == 0 and a == 1:
            del file
            file = np.loadtxt('./Data_dat/Diff_L_Photo.dat', unpack=True)
    if obj.name == "K+S0":
        plt.plot(x, y, label=f"$K^+\Sigma^0$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data/P4.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data/P5.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data/P6.dat', unpack=True)
        if Q2 == 0 and a == 1:
            del file
            file = np.loadtxt('./Data/Diff_S_Photo.dat', unpack=True)

    x_data = []
    y_data = []
    dy_data = []

    if a == 1 or b == 1 or len(file1) + len(file2) + len(file3) == 0:
        data_presence = False

    if data_presence:
        if len(file1) > 0:
            for i in range(0, len(file1[1])):
                if file1[1][i] == Q2 and file1[3][i] == cos_th:
                    x_data.append(file1[0][i])
                    y_data.append(file1[j][i])
                    dy_data.append(file1[j+1][i])
        if len(file2) > 0:
            for i in range(0, len(file2[1])):
                if file2[1][i] == Q2 and file2[3][i] == cos_th:
                    x_data.append(file2[0][i])
                    y_data.append(file2[j][i])
                    dy_data.append(file2[j+1][i])
        if len(file3) > 0:
            for i in range(0, len(file3[1])):
                if file3[1][i] == Q2 and file3[3][i] == cos_th:
                    x_data.append(file3[0][i])
                    y_data.append(file3[j][i])
                    dy_data.append(file3[j+1][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS data")

    if a == 1 and Q2 == 0:
        for i in range(0, len(file[1])):
            if abs(file[0][i] - W) < 1e-9:
                x_data.append(file[1][i])
                y_data.append(file[2][i])
                dy_data.append(file[3][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS Photo data")

    plt.xlabel('$\cos$')
    if name == "St":
        plt.ylabel('$S_{t}$, [nb/sr]')
        plt.title('$S_{t}$, [nb/sr]')
    if name == "Sl":
        plt.ylabel('$S_{l}$, [nb/sr]')
        plt.title('$S_{l}$, [nb/sr]')
    if name == "Slt":
        plt.ylabel('$S_{lt}$, [nb/sr]')
        plt.title('$S_{lt}$, [nb/sr]')
    if name == "Stt":
        plt.ylabel('$S_{tt}$, [nb/sr]')
        plt.title('$S_{tt}$, [nb/sr]')
    if name == "Su":
        plt.ylabel('$S_{u}$, [nb/sr]')
        plt.title('$S_{u}$, [nb/sr]')

    plt.legend(prop={'size': 20})

    # plt.show()

    if name == "St":
        plt.savefig(
            f"./pictures/St_cos_{obj.name}_{Q2}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Sl":
        plt.savefig(
            f"./pictures/Sl_cos_{obj.name}_{Q2}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Slt":
        plt.savefig(
            f"./pictures/Slt_cos_{obj.name}_{Q2}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Stt":
        plt.savefig(
            f"./pictures/Stt_cos_{obj.name}_{Q2}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Su":
        plt.savefig(
            f"./pictures/Su_cos_{obj.name}_{Q2}_{W}_{E_beam}.pdf", bbox_inches='tight')

    plt.clf()


def plot_str_Q2(name: str, obj: Model, W: float, cos_th: float, E_beam: float):

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0

    if name == "St":
        a = 1
    if name == "Sl":
        b = 1
    if name == "Slt":
        c = 1
        j = 6
    if name == "Stt":
        d = 1
        j = 8
    if name == "Su":
        e = 1
        j = 4

    obj.cos_th = cos_th
    obj.W = W
    obj.E_beam = E_beam

    x = np.arange(0.0, 5.0, 0.05)
    y = []
    dy = []

    for Q2 in x:
        obj.Q2 = Q2
        res = obj.Str_func_all()
        y.append(a*res[0] + b*res[2] + c*res[4] + d*res[6] + e*(res[0] + obj.eps(W, Q2)*res[2]))
        dy.append(a*res[1] + b*res[3] + c*res[5] + d*res[7] + e *
                  sqrt(res[1]**2 + obj.eps(W, Q2)**2*res[3]**2))

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    fig.set_dpi(100)

    plt.style.use('bmh')

    plt.scatter(x, y)
    plt.errorbar(x, y, yerr=dy, fmt="o",
                 label=f'W = {W} GeV\n $\cos$'+r'$\theta$ = ' + f'{cos_th}\nE = {E_beam} GeV')

    data_presence = True

    file1 = []
    file2 = []
    file3 = []
    file = []

    if obj.name == "K+L":
        plt.plot(x, y, label=f"$K^+\Lambda$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data_dat/P1.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data_dat/P2.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data_dat/P3.dat', unpack=True)

        file = np.loadtxt('./Data_dat/Diff_L_Photo.dat', unpack=True)

    if obj.name == "K+S0":
        plt.plot(x, y, label=f"$K^+\Sigma^0$")
        if E_beam == 2.567 or e == 0:
            file1 = np.loadtxt('./Data/P4.dat', unpack=True)
        if E_beam == 4.056 or e == 0:
            file2 = np.loadtxt('./Data/P5.dat', unpack=True)
        if E_beam == 5.499 or e == 0:
            file3 = np.loadtxt('./Data/P6.dat', unpack=True)

        file = np.loadtxt('./Data/Diff_S_Photo.dat', unpack=True)

    x_data = []
    y_data = []
    dy_data = []

    if a == 1 or b == 1 or len(file1) + len(file2) + len(file3) == 0:
        data_presence = False

    if data_presence:
        if len(file1) > 0:
            for i in range(0, len(file1[1])):
                if abs(file1[0][i] - W) < 1e-9 and abs(file1[3][i] - cos_th) < 1e-9:
                    x_data.append(file1[1][i])
                    y_data.append(file1[j][i])
                    dy_data.append(file1[j+1][i])
        if len(file2) > 0:
            for i in range(0, len(file2[1])):
                if abs(file2[0][i] - W) < 1e-9 and abs(file2[3][i] - cos_th) < 1e-9:
                    x_data.append(file2[1][i])
                    y_data.append(file2[j][i])
                    dy_data.append(file2[j+1][i])
        if len(file3) > 0:
            for i in range(0, len(file3[1])):
                if abs(file3[0][i] - W) < 1e-9 and abs(file3[3][i] - cos_th) < 1e-9:
                    x_data.append(file3[1][i])
                    y_data.append(file3[j][i])
                    dy_data.append(file3[j+1][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS data")

    if a == 1:
        for i in range(0, len(file[1])):
            if abs(file[0][i] - W) < 1e-9 and abs(file[1][i] - cos_th) < 1e-9:
                x_data.append(0)
                y_data.append(file[2][i])
                dy_data.append(file[3][i])
        if len(x_data) != 0:
            plt.scatter(x_data, y_data)
            plt.errorbar(x_data, y_data, yerr=dy_data, fmt="o", label=f"CLAS Photo data")

    plt.xlabel('$Q^2, GeV^{2}$')
    if name == "St":
        plt.ylabel('$S_{t}$, [nb/sr]')
        plt.title('$S_{t}$, [nb/sr]')
    if name == "Sl":
        plt.ylabel('$S_{l}$, [nb/sr]')
        plt.title('$S_{l}$, [nb/sr]')
    if name == "Slt":
        plt.ylabel('$S_{lt}$, [nb/sr]')
        plt.title('$S_{lt}$, [nb/sr]')
    if name == "Stt":
        plt.ylabel('$S_{tt}$, [nb/sr]')
        plt.title('$S_{tt}$, [nb/sr]')
    if name == "Su":
        plt.ylabel('$S_{u}$, [nb/sr]')
        plt.title('$S_{u}$, [nb/sr]')

    plt.legend(prop={'size': 20})

    # plt.show()

    if name == "St":
        plt.savefig(
            f"./pictures/St_Q2_{obj.name}_{cos_th}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Sl":
        plt.savefig(
            f"./pictures/Sl_Q2_{obj.name}_{cos_th}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Slt":
        plt.savefig(
            f"./pictures/Slt_Q2_{obj.name}_{cos_th}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Stt":
        plt.savefig(
            f"./pictures/Stt_Q2_{obj.name}_{cos_th}_{W}_{E_beam}.pdf", bbox_inches='tight')
    if name == "Su":
        plt.savefig(
            f"./pictures/Su_Q2_{obj.name}_{cos_th}_{W}_{E_beam}.pdf", bbox_inches='tight')

    plt.clf()
