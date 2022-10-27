from scipy import interpolate
from scipy.optimize import curve_fit
from numpy import cos, pi, sqrt, floor
from files import Files
from bin import Bin
import pandas as pd
import numpy as np


class Model:

    def __init__(self, options, W: list, Q2: list, cos_th: list, phi: list):
        # Check the channel
        assert options[
            'channel'] == "K+L" or options['channel'] == "K+S0", f"Channel name: {options['channel']} is invalid"
        # Check the phase space parameters
        for i in Q2:
            assert i >= 0, f"Q2: {i} is lower than zero"
        for i in cos_th:
            assert i >= -1 and i <= 1, f"cos_th: {i} is out of [-1, 1] interval"
        assert options['E_beam'] >= 0, f"E_beam: {options['E_beam']} is lower than zero"
        # if options['channel'] == "K+L":
        #    for i in W:
        #        assert i >= 1.61, f"W: {i} is less than threshold for KLambda channel"
        # elif options['channel'] == "K+S0":
        #    for i in W:
        #        assert i >= 1.68632, f"W: {i} is less than threshold for KSigma0 channel"

        # Assignment
        self.__name = options['channel']
        self.__W = W
        self.__Q2 = Q2
        self.__cos_th = cos_th
        self.__phi = phi
        self.__E_beam = options['E_beam']
        self.__Data = Files(options['channel'])
        self.__ratio_str = options['ratio_str']
        self.__add_factor = options['add_factor']
        self.__W_sys = options['W_sys']
        self.__err_option = options['err_option']

        vol = []
        vol_sf = []

        for i in W:
            for j in Q2:
                for k in cos_th:
                    vol_sf.append([i, j, k])
                    for z in phi:
                        vol.append([i, j, k, z])

        self.__req_volume = pd.DataFrame(vol, columns=['W', 'Q2', 'cos_th', 'phi'])
        self.__req_volume['E'] = self.__E_beam
        self.__req_volume['eps'] = self.eps(self.__req_volume['W'], self.__req_volume['Q2'])
        self.__req_volume_sf = pd.DataFrame(vol_sf, columns=['W', 'Q2', 'cos_th'])

        del vol, vol_sf

    @property
    def name(self):
        return self.__name

    @property
    def W(self):
        return self.__W

    @property
    def Q2(self):
        return self.__Q2

    @property
    def cos_th(self):
        return self.__cos_th

    @property
    def phi(self):
        return self.__phi

    @property
    def E_beam(self):
        return self.__E_beam

    @property
    def Data(self):
        return self.__Data

    @property
    def ratio_str(self):
        return self.__ratio_str

    @property
    def add_factor(self):
        return self.__add_factor

    @property
    def W_sys(self):
        return self.__W_sys

    @property
    def err_option(self):
        return self.__err_option

    @property
    def req_volume(self):
        return self.__req_volume

    @property
    def req_volume_sf(self):
        return self.__req_volume_sf

    @W.setter
    def W(self, value):
        self.__W = value

    @Q2.setter
    def Q2(self, value):
        self.__Q2 = value

    @cos_th.setter
    def cos_th(self, value):
        self.__cos_th = value

    @phi.setter
    def phi(self, value):
        self.__phi = value

    @E_beam.setter
    def E_beam(self, value):
        self.__E_beam = value

    @ratio_str.setter
    def ratio_str(self, value):
        self.__ratio_str = value

    @add_factor.setter
    def add_factor(self, value):
        self.__add_factor = value

    @W_sys.setter
    def W_sys(self, value):
        self.__W_sys = value

    @err_option.setter
    def err_option(self, value):
        self.__err_option = value

    def __repr__(self):
        return f'Model("{self.name}", {[x for x in self.W]}, {[x for x in self.Q2]}, {[x for x in self.cos_th]}, {[x for x in self.phi]}, {self.E_beam}, {self.ratio_str}, {self.add_factor}, {self.W_sys}, {self.err_option})'

    def __str__(self):

        print(f" Channel: {self.name}")
        print(f"\n Phase space parameters:")
        print(f"\t W\t= {self.W} GeV")
        print(f"\t Q2\t= {self.Q2} GeV2")
        print(f"\t cos_th\t= {self.cos_th}")
        print(f"\t phi\t= {self.phi} degree")
        print(f"\t E_beam\t= {self.E_beam} GeV\n\n")

        print(f" {self.ratio_str} - Special ratio method")
        print(f" {self.add_factor} - Factorized S_LT and S_TT")
        print(f" {self.W_sys} - Additional sys. error for W extrapolation for S_T str. func")

        if self.err_option == 3:
            print(f"\n Quad. extrapolation of errors for cos_th axis\n")
        elif self.err_option == 2:
            print(f"\n Linear exrapolation of errors(up to 100%) for cos_th axis\n")
        elif self.err_option == 1:
            print(f"\n Constant extrapolation of errors for cos_th axis\n")

        # print(self.Data)

        return ' '

    def __interp_cub_PP(self, data_table):
        result = pd.DataFrame()

        data_table = data_table.copy()

        data_table['dcs_plus'] = data_table['cs'] + data_table['dcs']

        W = data_table.iloc[0]['W']

        x = np.array(data_table['cos_th'])

        cs = np.array(data_table['cs'])
        dcs_plus = np.array(data_table['dcs_plus'])

        f1 = interpolate.interp1d(x, cs, kind='cubic', copy=False)
        f1_plus = interpolate.interp1d(x, dcs_plus, kind='cubic', copy=False)

        e1 = interpolate.InterpolatedUnivariateSpline(x, cs, k=2)
        e1_plus = interpolate.InterpolatedUnivariateSpline(x, dcs_plus, k=2)

        for cos_th in self.cos_th:
            if cos_th >= x.min() and cos_th <= x.max():
                result = pd.concat([result, pd.DataFrame({'W': [W], 'cos_th': [cos_th], 'dSt': [f1_plus(
                    cos_th).astype(np.float) - f1(cos_th).astype(np.float)]})], ignore_index=True)
            else:
                result = pd.concat([result, pd.DataFrame({'W': [W], 'cos_th': [cos_th], 'dSt': [e1_plus(
                    cos_th).astype(np.float) - e1(cos_th).astype(np.float)]})], ignore_index=True)

        return result

    def __interp_cub_PPSigma(self, data_table):
        result = pd.DataFrame()

        data_table = data_table.copy()

        data_table['dSigma_plus'] = data_table['Sigma'] + data_table['dSigma']

        W = data_table.iloc[0]['W']

        x = np.array(data_table['cos_th'])

        cs = np.array(data_table['Sigma'])
        dcs_plus = np.array(data_table['dSigma_plus'])

        f1 = interpolate.interp1d(x, cs, kind='cubic', copy=False)
        f1_plus = interpolate.interp1d(x, dcs_plus, kind='cubic', copy=False)

        e1 = interpolate.InterpolatedUnivariateSpline(x, cs, k=2)
        e1_plus = interpolate.InterpolatedUnivariateSpline(x, dcs_plus, k=2)

        for cos_th in self.cos_th:
            if cos_th >= x.min() and cos_th <= x.max():
                result = pd.concat([result, pd.DataFrame({'W': [W], 'cos_th': [cos_th], 'dSigma': [f1_plus(
                    cos_th).astype(np.float) - f1(cos_th).astype(np.float)]})], ignore_index=True)
            else:
                result = pd.concat([result, pd.DataFrame({'W': [W], 'cos_th': [cos_th], 'dSigma': [e1_plus(
                    cos_th).astype(np.float) - e1(cos_th).astype(np.float)]})], ignore_index=True)

        return result

    def __interp_cub(self, data_table):
        result = pd.DataFrame()

        data_table = data_table.copy()

        data_table['dSt_plus'] = data_table['St'] + data_table['dSt']
        data_table['dSl_plus'] = data_table['Sl'] + data_table['dSl']
        data_table['dSlt_plus'] = data_table['Slt'] + data_table['dSlt']
        data_table['dStt_plus'] = data_table['Stt'] + data_table['dStt']

        W = data_table.iloc[0]['W']
        Q2 = data_table.iloc[0]['Q2']

        x = np.array(data_table['cos_th'])

        St = np.array(data_table['St'])
        Sl = np.array(data_table['Sl'])
        Slt = np.array(data_table['Slt'])
        Stt = np.array(data_table['Stt'])

        dSt_plus = np.array(data_table['dSt_plus'])
        dSl_plus = np.array(data_table['dSl_plus'])
        dSlt_plus = np.array(data_table['dSlt_plus'])
        dStt_plus = np.array(data_table['dStt_plus'])

        type = ''

        if len(x) < 4:
            type = 'quadratic'
        else:
            type = 'cubic'

        f1 = interpolate.interp1d(x, St, kind=type, copy=False)
        f2 = interpolate.interp1d(x, Sl, kind=type, copy=False)
        f3 = interpolate.interp1d(x, Slt, kind=type, copy=False)
        f4 = interpolate.interp1d(x, Stt, kind=type, copy=False)

        f1_plus = interpolate.interp1d(x, dSt_plus, kind=type, copy=False)
        f2_plus = interpolate.interp1d(x, dSl_plus, kind=type, copy=False)
        f3_plus = interpolate.interp1d(x, dSlt_plus, kind=type, copy=False)
        f4_plus = interpolate.interp1d(x, dStt_plus, kind=type, copy=False)

        e1 = interpolate.InterpolatedUnivariateSpline(x, St, k=2)
        e2 = interpolate.InterpolatedUnivariateSpline(x, Sl, k=2)
        e3 = interpolate.InterpolatedUnivariateSpline(x, Slt, k=2)
        e4 = interpolate.InterpolatedUnivariateSpline(x, Stt, k=2)

        e1_plus = interpolate.InterpolatedUnivariateSpline(x, dSt_plus, k=2)
        e2_plus = interpolate.InterpolatedUnivariateSpline(x, dSl_plus, k=2)
        e3_plus = interpolate.InterpolatedUnivariateSpline(x, dSlt_plus, k=2)
        e4_plus = interpolate.InterpolatedUnivariateSpline(x, dStt_plus, k=2)

        for cos_th in self.cos_th:
            if cos_th >= x.min() and cos_th <= x.max():
                result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'dSt': [f1_plus(cos_th).astype(np.float) - f1(cos_th).astype(np.float)], 'dSl': [f2_plus(
                    cos_th).astype(np.float) - f2(cos_th).astype(np.float)], 'dSlt': [f3_plus(cos_th).astype(np.float) - f3(cos_th).astype(np.float)], 'dStt': [f4_plus(cos_th).astype(np.float) - f4(cos_th).astype(np.float)]})], ignore_index=True)
            else:
                result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'dSt': [e1_plus(cos_th).astype(np.float) - e1(cos_th).astype(np.float)], 'dSl': [e2_plus(
                    cos_th).astype(np.float) - e2(cos_th).astype(np.float)], 'dSlt': [e3_plus(cos_th).astype(np.float) - e3(cos_th).astype(np.float)], 'dStt': [e4_plus(cos_th).astype(np.float) - e4(cos_th).astype(np.float)]})], ignore_index=True)

        return result

    def __interp_cub_extraQ2(self, data_table, Q2_list):
        result = pd.DataFrame()

        data_table = data_table.copy()

        data_table['dSt_plus'] = data_table['St'] + data_table['dSt']
        data_table['dSl_plus'] = data_table['Sl'] + data_table['dSl']
        data_table['dSlt_plus'] = data_table['Slt'] + data_table['dSlt']
        data_table['dStt_plus'] = data_table['Stt'] + data_table['dStt']

        W = data_table.iloc[0]['W']
        cos_th = data_table.iloc[0]['cos_th']

        x = np.array(data_table['Q2'])

        St = np.array(data_table['St'])
        Sl = np.array(data_table['Sl'])
        Slt = np.array(data_table['Slt'])
        Stt = np.array(data_table['Stt'])

        dSt_plus = np.array(data_table['dSt_plus'])
        dSl_plus = np.array(data_table['dSl_plus'])
        dSlt_plus = np.array(data_table['dSlt_plus'])
        dStt_plus = np.array(data_table['dStt_plus'])

        e1 = interpolate.InterpolatedUnivariateSpline(x, St, k=2)
        e2 = interpolate.InterpolatedUnivariateSpline(x, Sl, k=2)
        e3 = interpolate.InterpolatedUnivariateSpline(x, Slt, k=2)
        e4 = interpolate.InterpolatedUnivariateSpline(x, Stt, k=2)

        e1_plus = interpolate.InterpolatedUnivariateSpline(x, dSt_plus, k=2)
        e2_plus = interpolate.InterpolatedUnivariateSpline(x, dSl_plus, k=2)
        e3_plus = interpolate.InterpolatedUnivariateSpline(x, dSlt_plus, k=2)
        e4_plus = interpolate.InterpolatedUnivariateSpline(x, dStt_plus, k=2)

        for Q2 in Q2_list:
            result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'dSt': [e1_plus(Q2).astype(np.float) - e1(Q2).astype(np.float)], 'dSl': [e2_plus(
                Q2).astype(np.float) - e2(Q2).astype(np.float)], 'dSlt': [e3_plus(Q2).astype(np.float) - e3(Q2).astype(np.float)], 'dStt': [e4_plus(Q2).astype(np.float) - e4(Q2).astype(np.float)]})], ignore_index=True)

        return result

    @staticmethod
    def __interp_cub_TT_LT(data_table, arg_th: float):

        x = np.array(data_table['cos_th'])
        if 'Slt' in data_table.columns:
            y = np.array(data_table['Slt'])
        elif 'Stt' in data_table.columns:
            y = np.array(data_table['Stt'])
        else:
            return 0

        x = np.append(x, -1)
        y = np.append(y, 0)

        f = interpolate.interp1d(x, y, kind='cubic', copy=False)

        del x, y

        return f(arg_th).astype(np.float)

    def __approx_cos_leg(self, data_table):
        W = data_table.iloc[0]['W']
        Q2 = data_table.iloc[0]['Q2']

        x = data_table['cos_th']

        St = data_table['St']
        dSt = data_table['dSt']
        Sl = data_table['Sl']
        dSl = data_table['dSl']
        Slt = data_table['Slt']
        dSlt = data_table['dSlt']
        Stt = data_table['Stt']
        dStt = data_table['dStt']

        def Su(x, a, b, c, d, e):
            y = a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8
            return y

        def SLT_f(x, a, b, c, d, e):
            y = np.sin(np.arccos(x))*(a + b*x + c*0.5*(3*x**2 - 1) + d *
                                      0.5*(5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8)
            return y

        def STT_f(x, a, b, c, d, e):
            y = (1-x**2)*(a + b*x + c*0.5*(3*x**2 - 1) + d*0.5 *
                          (5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8)
            return y

        def Su4(x, a, b, c, d):
            y = a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x)
            return y

        def SLT_f4(x, a, b, c, d):
            y = np.sin(np.arccos(x))*(a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x))
            return y

        def STT_f4(x, a, b, c, d):
            y = (1-x**2)*(a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x))
            return y

        def Su3(x, a, b, c):
            y = a + b*x + c*0.5*(3*x**2 - 1)
            return y

        def SLT_f3(x, a, b, c):
            y = np.sin(np.arccos(x))*(a + b*x + c*0.5*(3*x**2 - 1))
            return y

        def STT_f3(x, a, b, c):
            y = (1-x**2)*(a + b*x + c*0.5*(3*x**2 - 1))
            return y

        result = pd.DataFrame()

        if len(x) > 4:
            popt_St, pcov_St = curve_fit(Su, x, St, sigma=dSt, absolute_sigma=True)
            popt_Sl, pcov_Sl = curve_fit(Su, x, Sl, sigma=dSl, absolute_sigma=True)
            if self.add_factor:
                popt_Slt, pcov_Slt = curve_fit(SLT_f, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(STT_f, x, Stt, sigma=dStt, absolute_sigma=True)
            else:
                popt_Slt, pcov_Slt = curve_fit(Su, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(Su, x, Stt, sigma=dStt, absolute_sigma=True)

            for cos_th in self.cos_th:
                buff = {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th]}
                arg = cos_th

                buff['St'] = [Su(cos_th, popt_St[0], popt_St[1],
                                 popt_St[2], popt_St[3], popt_St[4])]
                buff['Sl'] = [Su(cos_th, popt_Sl[0], popt_Sl[1],
                                 popt_Sl[2], popt_Sl[3], popt_Sl[4])]

                if self.add_factor:
                    buff['Slt'] = [SLT_f(cos_th, popt_Slt[0], popt_Slt[1],
                                         popt_Slt[2], popt_Slt[3], popt_Slt[4])]
                    buff['Stt'] = [STT_f(cos_th, popt_Stt[0], popt_Stt[1],
                                         popt_Stt[2], popt_Stt[3], popt_Stt[4])]
                    if arg < x.min():
                        buff['Slt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Slt', 'dSlt']], arg)]
                        buff['Stt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Stt', 'dStt']], arg)]
                else:
                    buff['Slt'] = [Su(cos_th, popt_Slt[0], popt_Slt[1],
                                      popt_Slt[2], popt_Slt[3], popt_Slt[4])]
                    buff['Stt'] = [Su(cos_th, popt_Stt[0], popt_Stt[1],
                                      popt_Stt[2], popt_Stt[3], popt_Stt[4])]

                result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

            result_for_merge = self.__interp_cub(data_table)

            result_for_merge.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)

            result = result_for_merge.merge(result, on=['W', 'Q2', 'cos_th'])
            result.reset_index(inplace=True)

            if self.err_option == 1:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt']
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl']
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt']
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSt'] = data_table.iloc[len(data_table)-1]['dSt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSl'] = data_table.iloc[len(data_table)-1]['dSl']
                result.loc[result['cos_th'] >
                           x.max(), 'dSlt'] = data_table.iloc[len(data_table)-1]['dSlt']
                result.loc[result['cos_th'] >
                           x.max(), 'dStt'] = data_table.iloc[len(data_table)-1]['dStt']
            if self.err_option == 2:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)

                result.loc[result['cos_th'] > x.max(), 'dSt'] = data_table.iloc[len(
                    data_table)-1]['dSt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSl'] = data_table.iloc[len(
                    data_table)-1]['dSl']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
            result.loc[result['cos_th'] < x.min(), 'St'] = Su(
                x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
            result.loc[result['cos_th'] < x.min(), 'Sl'] = Su(
                x.min(), popt_Sl[0], popt_Sl[1], popt_Sl[2], popt_Sl[3], popt_Sl[4])
            result.loc[result['cos_th'] > x.max(), 'St'] = Su(
                x.max(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
            result.loc[result['cos_th'] > x.max(), 'Sl'] = Su(
                x.max(), popt_Sl[0], popt_Sl[1], popt_Sl[2], popt_Sl[3], popt_Sl[4])
            if not(self.add_factor):
                result.loc[result['cos_th'] < x.min(), 'Slt'] = Su(
                    x.min(), popt_Slt[0], popt_Slt[1], popt_Slt[2], popt_Slt[3], popt_Slt[4])
                result.loc[result['cos_th'] < x.min(), 'Stt'] = Su(
                    x.min(), popt_Stt[0], popt_Stt[1], popt_Stt[2], popt_Stt[3], popt_Stt[4])
                result.loc[result['cos_th'] > x.max(), 'Slt'] = Su(
                    x.max(), popt_Slt[0], popt_Slt[1], popt_Slt[2], popt_Slt[3], popt_Slt[4])
                result.loc[result['cos_th'] > x.max(), 'Stt'] = Su(
                    x.max(), popt_Stt[0], popt_Stt[1], popt_Stt[2], popt_Stt[3], popt_Stt[4])
            if self.add_factor:
                result.loc[(result['cos_th'] < x.min()), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] < x.min()), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] > x.max()), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())
                result.loc[(result['cos_th'] > x.max()), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())

            result.loc[result['St'] < 0, 'St'] = 0
            result.loc[result['Sl'] < 0, 'Sl'] = 0

        elif len(x) == 4:
            popt_St, pcov_St = curve_fit(Su4, x, St, sigma=dSt, absolute_sigma=True)
            popt_Sl, pcov_Sl = curve_fit(Su4, x, Sl, sigma=dSl, absolute_sigma=True)
            if self.add_factor:
                popt_Slt, pcov_Slt = curve_fit(SLT_f4, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(STT_f4, x, Stt, sigma=dStt, absolute_sigma=True)
            else:
                popt_Slt, pcov_Slt = curve_fit(Su4, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(Su4, x, Stt, sigma=dStt, absolute_sigma=True)

            for cos_th in self.cos_th:
                buff = {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th]}
                arg = cos_th

                buff['St'] = [Su4(cos_th, popt_St[0], popt_St[1], popt_St[2], popt_St[3])]
                buff['Sl'] = [Su4(cos_th, popt_Sl[0], popt_Sl[1], popt_Sl[2], popt_Sl[3])]

                if self.add_factor:
                    buff['Slt'] = [SLT_f4(cos_th, popt_Slt[0], popt_Slt[1],
                                          popt_Slt[2], popt_Slt[3])]
                    buff['Stt'] = [STT_f4(cos_th, popt_Stt[0], popt_Stt[1],
                                          popt_Stt[2], popt_Stt[3])]
                    if arg < x.min():
                        buff['Slt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Slt', 'dSlt']], arg)]
                        buff['Stt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Stt', 'dStt']], arg)]
                else:
                    buff['Slt'] = [Su4(cos_th, popt_Slt[0], popt_Slt[1], popt_Slt[2], popt_Slt[3])]
                    buff['Stt'] = [Su4(cos_th, popt_Stt[0], popt_Stt[1], popt_Stt[2], popt_Stt[3])]

                result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

            result_for_merge = self.__interp_cub(data_table)

            result_for_merge.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)

            result = result_for_merge.merge(result, on=['W', 'Q2', 'cos_th'])
            result.reset_index(inplace=True)

            if self.err_option == 1:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt']
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl']
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt']
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSt'] = data_table.iloc[len(data_table)-1]['dSt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSl'] = data_table.iloc[len(data_table)-1]['dSl']
                result.loc[result['cos_th'] >
                           x.max(), 'dSlt'] = data_table.iloc[len(data_table)-1]['dSlt']
                result.loc[result['cos_th'] >
                           x.max(), 'dStt'] = data_table.iloc[len(data_table)-1]['dStt']
            if self.err_option == 2:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)

                result.loc[result['cos_th'] > x.max(), 'dSt'] = data_table.iloc[len(
                    data_table)-1]['dSt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSl'] = data_table.iloc[len(
                    data_table)-1]['dSl']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
            result.loc[result['cos_th'] < x.min(), 'St'] = Su4(
                x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3])
            result.loc[result['cos_th'] < x.min(), 'Sl'] = Su4(
                x.min(), popt_Sl[0], popt_Sl[1], popt_Sl[2], popt_Sl[3])
            result.loc[result['cos_th'] > x.max(), 'St'] = Su4(
                x.max(), popt_St[0], popt_St[1], popt_St[2], popt_St[3])
            result.loc[result['cos_th'] > x.max(), 'Sl'] = Su4(
                x.max(), popt_Sl[0], popt_Sl[1], popt_Sl[2], popt_Sl[3])
            if not(self.add_factor):
                result.loc[result['cos_th'] < x.min(), 'Slt'] = Su4(
                    x.min(), popt_Slt[0], popt_Slt[1], popt_Slt[2], popt_Slt[3])
                result.loc[result['cos_th'] < x.min(), 'Stt'] = Su4(
                    x.min(), popt_Stt[0], popt_Stt[1], popt_Stt[2], popt_Stt[3])
                result.loc[result['cos_th'] > x.max(), 'Slt'] = Su4(
                    x.max(), popt_Slt[0], popt_Slt[1], popt_Slt[2], popt_Slt[3])
                result.loc[result['cos_th'] > x.max(), 'Stt'] = Su4(
                    x.max(), popt_Stt[0], popt_Stt[1], popt_Stt[2], popt_Stt[3])
            if self.add_factor:
                result.loc[(result['cos_th'] < x.min()), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] < x.min()), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] > x.max()), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())
                result.loc[(result['cos_th'] > x.max()), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())

            result.loc[result['St'] < 0, 'St'] = 0
            result.loc[result['Sl'] < 0, 'Sl'] = 0

        elif len(x) < 4:
            popt_St, pcov_St = curve_fit(Su3, x, St, sigma=dSt, absolute_sigma=True)
            popt_Sl, pcov_Sl = curve_fit(Su3, x, Sl, sigma=dSl, absolute_sigma=True)
            if self.add_factor:
                popt_Slt, pcov_Slt = curve_fit(SLT_f3, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(STT_f3, x, Stt, sigma=dStt, absolute_sigma=True)
            else:
                popt_Slt, pcov_Slt = curve_fit(Su3, x, Slt, sigma=dSlt, absolute_sigma=True)
                popt_Stt, pcov_Stt = curve_fit(Su3, x, Stt, sigma=dStt, absolute_sigma=True)

            for cos_th in self.cos_th:
                buff = {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th]}
                arg = cos_th

                buff['St'] = [Su3(cos_th, popt_St[0], popt_St[1], popt_St[2])]
                buff['Sl'] = [Su3(cos_th, popt_Sl[0], popt_Sl[1], popt_Sl[2])]

                if self.add_factor:
                    buff['Slt'] = [SLT_f3(cos_th, popt_Slt[0], popt_Slt[1], popt_Slt[2])]
                    buff['Stt'] = [STT_f3(cos_th, popt_Stt[0], popt_Stt[1], popt_Stt[2])]
                    if arg < x.min():
                        buff['Slt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Slt', 'dSlt']], arg)]
                        buff['Stt'] = [self.__interp_cub_TT_LT(
                            data_table[['cos_th', 'Stt', 'dStt']], arg)]
                else:
                    buff['Slt'] = [Su3(cos_th, popt_Slt[0], popt_Slt[1], popt_Slt[2])]
                    buff['Stt'] = [Su3(cos_th, popt_Stt[0], popt_Stt[1], popt_Stt[2])]

                result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

            result_for_merge = self.__interp_cub(data_table)

            result_for_merge.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)

            result = result_for_merge.merge(result, on=['W', 'Q2', 'cos_th'])
            result.reset_index(inplace=True)

            if self.err_option == 1:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt']
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl']
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt']
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSt'] = data_table.iloc[len(data_table)-1]['dSt']
                result.loc[result['cos_th'] >
                           x.max(), 'dSl'] = data_table.iloc[len(data_table)-1]['dSl']
                result.loc[result['cos_th'] >
                           x.max(), 'dSlt'] = data_table.iloc[len(data_table)-1]['dSlt']
                result.loc[result['cos_th'] >
                           x.max(), 'dStt'] = data_table.iloc[len(data_table)-1]['dStt']
            if self.err_option == 2:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = data_table.iloc[0]['dSt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSl'] = data_table.iloc[0]['dSl'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] < x.min(), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)

                result.loc[result['cos_th'] > x.max(), 'dSt'] = data_table.iloc[len(
                    data_table)-1]['dSt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSl'] = data_table.iloc[len(
                    data_table)-1]['dSl']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)
            result.loc[result['cos_th'] < x.min(), 'St'] = Su3(
                x.min(), popt_St[0], popt_St[1], popt_St[2])
            result.loc[result['cos_th'] < x.min(), 'Sl'] = Su3(
                x.min(), popt_Sl[0], popt_Sl[1], popt_Sl[2])
            result.loc[result['cos_th'] > x.max(), 'St'] = Su3(
                x.max(), popt_St[0], popt_St[1], popt_St[2])
            result.loc[result['cos_th'] > x.max(), 'Sl'] = Su3(
                x.max(), popt_Sl[0], popt_Sl[1], popt_Sl[2])
            if not(self.add_factor):
                result.loc[result['cos_th'] < x.min(), 'Slt'] = Su3(
                    x.min(), popt_Slt[0], popt_Slt[1], popt_Slt[2])
                result.loc[result['cos_th'] < x.min(), 'Stt'] = Su3(
                    x.min(), popt_Stt[0], popt_Stt[1], popt_Stt[2])
                result.loc[result['cos_th'] > x.max(), 'Slt'] = Su3(
                    x.max(), popt_Slt[0], popt_Slt[1], popt_Slt[2])
                result.loc[result['cos_th'] > x.max(), 'Stt'] = Su3(
                    x.max(), popt_Stt[0], popt_Stt[1], popt_Stt[2])
            if self.add_factor:
                result.loc[(result['cos_th'] < x.min()), 'dSlt'] = data_table.iloc[0]['dSlt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] < x.min()), 'dStt'] = data_table.iloc[0]['dStt'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] > x.max()), 'dSlt'] = data_table.iloc[len(
                    data_table)-1]['dSlt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())
                result.loc[(result['cos_th'] > x.max()), 'dStt'] = data_table.iloc[len(
                    data_table)-1]['dStt']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())

            result.loc[result['St'] < 0, 'St'] = 0
            result.loc[result['Sl'] < 0, 'Sl'] = 0

        del x, St, dSt, Sl, dSl, Slt, dSlt, Stt, dStt
        del buff, popt_St, pcov_St, popt_Sl, pcov_Sl
        del popt_Slt, pcov_Slt, popt_Stt, pcov_Stt
        del result_for_merge

        return result

    def eps(self, W: float, Q2: float):
        m_p = 0.93827
        nu = (W*W + Q2 - m_p*m_p)/(2*m_p)
        return 1/(1 + 2*(nu*nu + Q2)/(4*(self.E_beam - nu)*self.E_beam - Q2))

    @staticmethod
    def eps_beam(W: float, Q2: float, E0: float):
        m_p = 0.93827
        nu = (W*W + Q2 - m_p*m_p)/(2*m_p)
        return 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2))

    def __isData1(self, W: float, Q2: float):
        if self.name == "K+L":
            return Q2 <= 1.0 and Q2 >= 0.65 and W <= 1.975 and W >= 1.65
        elif self.name == "K+S0":
            return Q2 <= 1.0 and Q2 >= 0.65 and W <= 2.05 and W >= 1.725

        return False

    def __isData2(self, W: float, Q2: float):
        if self.name == "K+L":
            return (Q2 <= 2.55 and Q2 >= 1.0 and W <= 2.15 and W >= 1.65) or (Q2 <= 2.05 and Q2 >= 1.0 and W <= 2.35 and W >= 2.15)
        elif self.name == "K+S0":
            return (Q2 <= 2.55 and Q2 >= 1.0 and W <= 2.15 and W >= 1.75) or (Q2 <= 2.05 and Q2 >= 1.0 and W <= 2.35 and W >= 2.15)

        return False

    def __isData3(self, W: float, Q2: float):
        if self.name == "K+L":
            return Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.575 and W >= 1.63
        elif self.name == "K+S0":
            return Q2 <= 3.45 and Q2 >= 1.8 and W <= 2.575 and W >= 1.695

        return False

    def __isW1(self, W: float):
        if self.name == "K+L":
            return W <= 1.975 and W >= 1.65
        elif self.name == "K+S0":
            return W <= 2.05 and W >= 1.725

        return False

    def __isW2(self, W: float):
        if self.name == "K+L":
            return W >= 1.65 and W <= 2.35
        elif self.name == "K+S0":
            return W >= 1.75 and W <= 2.35

        return False

    def __isW3(self, W: float):
        if self.name == "K+L":
            return W <= 2.575 and W >= 1.63
        elif self.name == "K+S0":
            return W <= 2.575 and W >= 1.695

        return False

    def __isWforQ2inter(self, W: float, Q2: float):
        if Q2 < 1.8:
            if self.name == "K+L":
                return W <= 2.575 and W >= 1.65
            else:
                if Q2 > 1.0:
                    return W <= 2.575 and W >= 1.75
                else:
                    return W <= 2.575 and W >= 1.725
        else:
            if self.name == "K+L":
                return W <= 2.575 and W >= 1.63
            else:
                return W <= 2.575 and W >= 1.695

    def __isSigma(self, W: float):
        if self.name == "K+L":
            return W <= 2.18 and W >= 1.72
        else:
            return W <= 2.17 and W >= 1.78

    def __Give_MinMaxE(self, W: float, Q2: float, zone: int):
        result = {}

        result["W"] = [W]
        result["Q2"] = [Q2]

        if zone == 1:
            result["W_min"] = float(floor((W - 0.025)*20))/20 + 0.025
            result["W_max"] = result["W_min"] + 0.05

            if W < 1.725:
                result["W_min"] = 1.65
                result["W_max"] = 1.725

            if W > 1.925:
                result["W_min"] = 1.925
                result["W_max"] = 1.975

            result["Q2_min"] = 0.65
            result["Q2_max"] = 1.0

            result["E"] = 2.567
        elif zone == 2:
            result["Q2_min"] = float(floor((Q2 - 0.05)*2))/2 + 0.05
            result["Q2_max"] = result["Q2_min"] + 0.5

            if Q2 < 1.55:
                result["Q2_min"] = 1.0
                result["Q2_max"] = 1.55

            if Q2 > 2.05:
                result["Q2_min"] = 2.05
                result["Q2_max"] = 2.55

            result["W_min"] = int(1000*(float(floor((W - 0.05)*10))/10 + 0.05))/1000.0
            result["W_max"] = result["W_min"] + 0.1

            if abs(W - 2.15) < 1e-9:
                result["W_min"] = 2.05
                result["W_max"] = 2.15

            if abs(W - 1.65) < 1e-9:
                result["W_min"] = 1.65
                result["W_max"] = 1.75

            if abs(W - 1.75) < 1e-9:
                result["W_min"] = 1.75
                result["W_max"] = 1.85

            if abs(W - 2.35) < 1e-9:
                result["W_min"] = 2.25
                result["W_max"] = 2.35

            result["E"] = 4.056
        elif zone == 3:
            result["W_min"] = float(floor((W - 0.025)*20))/20 + 0.025
            result["W_max"] = result["W_min"] + 0.05

            if W < 1.675 and self.name == "K+L":
                result["W_min"] = 1.63
                result["W_max"] = 1.675

            if W < 1.725 and self.name == "K+S0":
                result["W_min"] = 1.695
                result["W_max"] = 1.725

            if W > 2.525:
                result["W_min"] = 2.525
                result["W_max"] = 2.575

            if Q2 < 2.6:
                result["Q2_min"] = 1.8
                result["Q2_max"] = 2.6

            if 2.6 <= Q2:
                result["Q2_min"] = 2.6
                result["Q2_max"] = 3.45

            result["E"] = 5.499

        result["W_min"] = float("{0:.3f}".format(result["W_min"]))
        result["W_max"] = float("{0:.3f}".format(result["W_max"]))
        result["Q2_min"] = float("{0:.2f}".format(result["Q2_min"]))
        result["Q2_max"] = float("{0:.2f}".format(result["Q2_max"]))

        W_min = result["W_min"]
        W_max = result["W_max"]
        Q2_min = result["Q2_min"]
        Q2_max = result["Q2_max"]

        result = [{'W0': W, 'Q20': Q2, 'W': result["W_min"], 'Q2': result["Q2_min"], 'E': result["E"]}, {'W0': W, 'Q20': Q2, 'W': result["W_max"], 'Q2': result["Q2_min"], 'E': result["E"]}, {
            'W0': W, 'Q20': Q2, 'W': result["W_min"], 'Q2': result["Q2_max"], 'E': result["E"]}, {'W0': W, 'Q20': Q2, 'W': result["W_max"], 'Q2': result["Q2_max"], 'E': result["E"]}]

        W_Q2 = pd.DataFrame({'W0': [W], 'Q20': [Q2], 'W_min': [W_min],
                            'Q2_min': [Q2_min], 'W_max': [W_max], 'Q2_max': [Q2_max]})

        return result, W_Q2

    def __giveData(self, data_table, zone: int):

        buff = pd.DataFrame()
        W_Q2 = pd.DataFrame()

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            prior, minor = self.__Give_MinMaxE(W, Q2, zone)
            buff = pd.concat([buff, pd.DataFrame(prior)])
            W_Q2 = pd.concat([W_Q2, minor], ignore_index=True)

        result = buff.drop(columns=['E'])
        buff = buff.drop(columns=['W0', 'Q20'])
        buff = buff.groupby(['W', 'Q2'], as_index=False).mean()
        out = pd.DataFrame()

        for index, row in buff.iterrows():
            out = pd.concat([out, self.__approx_cos_leg(self.Data.Data.loc[(abs(self.Data.Data['W'] - row['W']) < 1e-6) & (
                abs(self.Data.Data['Q2'] - row['Q2']) < 1e-6) & (abs(self.Data.Data['E'] - row['E']) < 1e-5)])])

        result.set_index(['W', 'Q2'], inplace=True, drop=True)
        out.set_index(['W', 'Q2'], inplace=True, drop=True)

        result = out.merge(result, on=['W', 'Q2'])  # Check all indecies
        result.reset_index(inplace=True)

        result_final = pd.DataFrame()

        W_Q2.set_index(['W_min', 'Q2_min', 'W_max', 'Q2_max'], inplace=True, drop=True)

        W_Q2_dict = {}

        for index, row in W_Q2.iterrows():
            W_Q2_dict.setdefault(index, []).append([row['W0'], row['Q20']])

        for i in W_Q2_dict.keys():
            for cos_th in self.cos_th:
                del buff
                buff = pd.DataFrame()
                buff = pd.concat([buff, (result.loc[(result['W'] == i[0]) & (result['Q2'] == i[1]) & (
                    result['cos_th'] == cos_th)]).drop(columns=['W0', 'Q20'])], ignore_index=True)
                buff = pd.concat([buff, (result.loc[(result['W'] == i[0]) & (result['Q2'] == i[3]) & (
                    result['cos_th'] == cos_th)]).drop(columns=['W0', 'Q20'])], ignore_index=True)
                buff = pd.concat([buff, (result.loc[(result['W'] == i[2]) & (result['Q2'] == i[1]) & (
                    result['cos_th'] == cos_th)]).drop(columns=['W0', 'Q20'])], ignore_index=True)
                buff = pd.concat([buff, (result.loc[(result['W'] == i[2]) & (result['Q2'] == i[3]) & (
                    result['cos_th'] == cos_th)]).drop(columns=['W0', 'Q20'])], ignore_index=True)

                buff.drop_duplicates(subset=['W', 'Q2', 'cos_th'], inplace=True)
                buff.reset_index(inplace=True)

                x = np.array(buff['W'])
                y = np.array(buff['Q2'])

                f1 = np.array(buff['St'])
                df1 = np.array(buff['dSt'])
                f2 = np.array(buff['Sl'])
                df2 = np.array(buff['dSl'])
                f3 = np.array(buff['Slt'].astype(np.float))
                df3 = np.array(buff['dSlt'])
                f4 = np.array(buff['Stt'].astype(np.float))
                df4 = np.array(buff['dStt'])

                St_int = interpolate.interp2d(x, y, f1, kind='linear', copy=False)
                dSt_int = interpolate.interp2d(x, y, df1, kind='linear', copy=False)
                Sl_int = interpolate.interp2d(x, y, f2,  kind='linear', copy=False)
                dSl_int = interpolate.interp2d(x, y, df2, kind='linear', copy=False)
                Slt_int = interpolate.interp2d(x, y, f3,  kind='linear', copy=False)
                dSlt_int = interpolate.interp2d(x, y, df3,  kind='linear', copy=False)
                Stt_int = interpolate.interp2d(x, y, f4, kind='linear', copy=False)
                dStt_int = interpolate.interp2d(x, y, df4, kind='linear', copy=False)

                for [W, Q2] in W_Q2_dict[i]:
                    St = St_int(W, Q2).astype(np.float)[0]
                    dSt = dSt_int(W, Q2).astype(np.float)[0]
                    Sl = Sl_int(W, Q2).astype(np.float)[0]
                    dSl = dSl_int(W, Q2).astype(np.float)[0]
                    Slt = Slt_int(W, Q2).astype(np.float)[0]
                    dSlt = dSlt_int(W, Q2).astype(np.float)[0]
                    Stt = Stt_int(W, Q2).astype(np.float)[0]
                    dStt = dStt_int(W, Q2).astype(np.float)[0]

                    df_new_row = pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [St], 'dSt': [dSt], 'Sl': [
                                              Sl], 'dSl': [dSl], 'Slt': [Slt], 'dSlt': [dSlt], 'Stt': [Stt], 'dStt': [dStt]})

                    result_final = pd.concat([result_final, df_new_row], ignore_index=True)

        result_final.loc[result_final['St'] < 0, 'St'] = 0
        result_final.loc[result_final['Sl'] < 0, 'Sl'] = 0

        del buff, result, out, df_new_row, W_Q2, W_Q2_dict

        return result_final

    def __Photo_diff_fixed(self, data_table):
        reader = pd.DataFrame()
        result_final = pd.DataFrame()
        result_for_merge = pd.DataFrame()

        buff = {}

        def Su(x, a, b, c, d, e):
            y = a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8
            return y

        for W in data_table:
            W = float("{0:.3f}".format(W))
            reader = self.Data.Data_Diff.loc[self.Data.Data_Diff['W'] == W]

            x = reader['cos_th']
            St = reader['cs']
            dSt = reader['dcs']

            popt_St, pcov_St = curve_fit(Su, x, St, sigma=dSt, absolute_sigma=True)

            result = pd.DataFrame()

            for cos_th in self.cos_th:
                buff = {'W': [W], 'cos_th': [cos_th]}
                buff['St'] = [Su(cos_th, popt_St[0], popt_St[1],
                                 popt_St[2], popt_St[3], popt_St[4])]

                result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

            result_for_merge = self.__interp_cub_PP(reader)

            result_for_merge.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result = result_for_merge.merge(result, on=['W', 'cos_th'])
            result.reset_index(inplace=True)

            if self.err_option == 1:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = reader.iloc[0]['dcs']
                result.loc[result['cos_th'] >
                           x.max(), 'dSt'] = reader.iloc[len(reader)-1]['dcs']

            if self.err_option == 2:
                result.loc[result['cos_th'] < x.min(), 'dSt'] = reader.iloc[0]['dcs'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSt'] = reader.iloc[len(
                    reader)-1]['dcs']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)

            result.loc[result['cos_th'] < x.min(), 'St'] = Su(
                x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])

            result.loc[result['St'] < 0, 'St'] = 0

            result_final = pd.concat([result_final, result], ignore_index=True)
            del result

        del reader, result_for_merge, buff

        return result_final

    def __Photo_diff(self, data_table):
        W_list = data_table['W']

        result, buff = pd.DataFrame(), pd.DataFrame()

        for W in W_list:
            W_min = float(int(floor(W*100)))/100 + 0.005

            if W < W_min:
                W_max = W_min
                W_min = W_max - 0.01
            else:
                W_max = W_min + 0.01

            if abs(W_min - 1.955) < 1e-9:
                W_min = 1.945
            if abs(W_max - 1.955) < 1e-9:
                W_max = 1.965

            if abs(W_min - 2.735) < 1e-9:
                W_min = 2.725
            if abs(W_max - 2.735) < 1e-9:
                W_max = 2.755

            if abs(W_min - 2.745) < 1e-9:
                W_min = 2.725
            if abs(W_max - 2.745) < 1e-9:
                W_max = 2.755

            W_min = float("{0:.3f}".format(W_min))
            W_max = float("{0:.3f}".format(W_max))

            result = pd.concat([pd.DataFrame({'W0': [W], 'W': [W_min]}), result], ignore_index=True)
            result = pd.concat([pd.DataFrame({'W0': [W], 'W': [W_max]}), result], ignore_index=True)

        buff = result['W']
        buff.drop_duplicates(inplace=True)

        result_for_merge = self.__Photo_diff_fixed(buff)

        result_for_merge.set_index('W', inplace=True, drop=True)
        result = result.set_index('W', drop=True)
        result = result.merge(result_for_merge, on=['W'])
        result.reset_index(inplace=True)

        result_final = pd.DataFrame()

        for W in W_list:
            for cos_th in self.cos_th:
                buff = (result.loc[(result['W0'] == W) & (result['cos_th'] == cos_th)]).drop(
                    columns=['W0', 'cos_th'])
                buff.reset_index(inplace=True, drop=True)

                x = np.array(buff['W'])

                f1 = np.array(buff['St'])
                df1 = np.array(buff['dSt'])

                St_int = interpolate.interp1d(x, f1, kind='linear', copy=False)
                dSt_int = interpolate.interp1d(x, df1, kind='linear', copy=False)

                St = St_int(W).astype(np.float)
                dSt = dSt_int(W).astype(np.float)

                result_final = pd.concat([result_final, pd.DataFrame(
                    {'W': [W], 'Q2': [0], 'cos_th': [cos_th], 'St': [St], 'dSt': [dSt]})], ignore_index=True)

        result_final.loc[result_final['St'] < 0, 'St'] = 0

        result_final['Sl'] = 0
        result_final['dSl'] = 0
        result_final['Slt'] = 0
        result_final['dSlt'] = 0

        del buff, W_list, result, result_for_merge

        return result_final

    def __Photo_Sigma_fixed(self, data_table):
        reader = pd.DataFrame()
        result_final = pd.DataFrame()
        result_for_merge = pd.DataFrame()

        buff = {}

        def Su(x, a, b, c, d, e):
            y = a + b*x + c*0.5*(3*x**2 - 1) + d*0.5*(5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8
            return y

        def STT_f(x, a, b, c, d, e):
            y = (1-x**2)*(a + b*x + c*0.5*(3*x**2 - 1) + d*0.5 *
                          (5*x**3 - 3*x) + e*(35*x**4 - 30*x**2 + 3)/8)
            return y

        for W in data_table:
            W = float("{0:.3f}".format(W))
            reader = self.Data.Data_Sigma.loc[self.Data.Data_Sigma['W'] == W]

            x = reader['cos_th']
            St = reader['Sigma']
            dSt = reader['dSigma']

            if self.add_factor:
                popt_St, pcov_St = curve_fit(STT_f, x, St, sigma=dSt, absolute_sigma=True)
            else:
                popt_St, pcov_St = curve_fit(Su, x, St, sigma=dSt, absolute_sigma=True)

            result = pd.DataFrame()

            for cos_th in self.cos_th:
                buff = {'W': [W], 'cos_th': [cos_th]}
                if self.add_factor:
                    buff['Sigma'] = [STT_f(cos_th, popt_St[0], popt_St[1],
                                           popt_St[2], popt_St[3], popt_St[4])]
                else:
                    buff['Sigma'] = [Su(cos_th, popt_St[0], popt_St[1],
                                        popt_St[2], popt_St[3], popt_St[4])]

                result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

            result_for_merge = self.__interp_cub_PPSigma(reader)

            result_for_merge.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result = result_for_merge.merge(result, on=['W', 'cos_th'])
            result.reset_index(inplace=True)

            if self.err_option == 1:
                result.loc[result['cos_th'] < x.min(), 'dSigma'] = reader.iloc[0]['dSigma']
                result.loc[result['cos_th'] >
                           x.max(), 'dSigma'] = reader.iloc[len(reader)-1]['dSigma']

            if self.err_option == 2:
                result.loc[result['cos_th'] < x.min(), 'dSigma'] = reader.iloc[0]['dSigma'] * \
                    ((-result.loc[(result['cos_th'] < x.min()),
                     'cos_th'] + x.min())/(1 + x.min()) + 1)
                result.loc[result['cos_th'] > x.max(), 'dSigma'] = reader.iloc[len(
                    reader)-1]['dSigma']*((result.loc[(result['cos_th'] > x.max()), 'cos_th'] - x.max())/(1 - x.max()) + 1)

            if not(self.add_factor):
                result.loc[result['cos_th'] < x.min(), 'Sigma'] = Su(
                    x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
                result.loc[result['cos_th'] < x.min(), 'Sigma'] = Su(
                    x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
                result.loc[result['cos_th'] > x.max(), 'Sigma'] = Su(
                    x.max(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
            if self.add_factor:
                result.loc[result['cos_th'] < x.min(), 'Sigma'] = STT_f(
                    x.min(), popt_St[0], popt_St[1], popt_St[2], popt_St[3], popt_St[4])
                result.loc[(result['cos_th'] < x.min()), 'dSigma'] = reader.iloc[0]['dSigma'] * \
                    (result.loc[(result['cos_th'] < x.min()), 'cos_th'] + 1)/(1 + x.min())
                result.loc[(result['cos_th'] > x.max()), 'dSigma'] = reader.iloc[len(
                    reader)-1]['dSigma']*(result.loc[(result['cos_th'] > x.max()), 'cos_th'] - 1)/(-1 + x.max())

            result.loc[result['Sigma'] < 0, 'Sigma'] = 0

            result_final = pd.concat([result_final, result], ignore_index=True)
            del result

        del reader, result_for_merge, buff

        return result_final

    def __Photo_Sigma(self, data_table):
        W_list = data_table['W']

        result, buff = pd.DataFrame(), pd.DataFrame()

        for W in W_list:
            if self.name == "K+L":
                W_min = float(floor(W*50))/50
                W_max = W_min + 0.02
                if W_max > 2.18:
                    W_max = 2.18
                    W_min = 2.16
            elif self.name == "K+S0":
                W_min = float(floor((W-0.01)*25))/25 + 0.01
                W_max = W_min + 0.04

                if W >= 1.78 and W < 1.81:
                    W_min = 1.78
                    W_max = 1.81

                if W > 2.13:
                    W_min = 2.13
                    W_max = 2.17

            W_min = float("{0:.3f}".format(W_min))
            W_max = float("{0:.3f}".format(W_max))

            result = pd.concat(
                [pd.DataFrame({'W0': [W], 'W': [W_min]}), result], ignore_index=True)
            result = pd.concat(
                [pd.DataFrame({'W0': [W], 'W': [W_max]}), result], ignore_index=True)

        buff = result['W']
        buff.drop_duplicates(inplace=True)

        result_for_merge = self.__Photo_Sigma_fixed(buff)

        result_for_merge.set_index('W', inplace=True, drop=True)
        result = result.set_index('W', drop=True)
        result = result.merge(result_for_merge, on=['W'])
        result.reset_index(inplace=True)

        result_final = pd.DataFrame()

        for W in W_list:
            for cos_th in self.cos_th:
                buff = (result.loc[(result['W0'] == W) & (result['cos_th'] == cos_th)]).drop(
                    columns=['W0', 'cos_th'])
                buff.reset_index(inplace=True, drop=True)

                x = np.array(buff['W'])

                f1 = np.array(buff['Sigma'])
                df1 = np.array(buff['dSigma'])

                St_int = interpolate.interp1d(x, f1, kind='linear', copy=False)
                dSt_int = interpolate.interp1d(x, df1, kind='linear', copy=False)

                St = St_int(W).astype(np.float)
                dSt = dSt_int(W).astype(np.float)

                result_final = pd.concat([result_final, pd.DataFrame(
                    {'W': [W], 'Q2': [0], 'cos_th': [cos_th], 'Sigma': [St], 'dSigma': [dSt]})], ignore_index=True)

        result_final.loc[result_final['Sigma'] < 0, 'Sigma'] = 0

        del buff, W_list, result, result_for_merge

        return result_final

    def __Str_func_for_lower(self, data_table):
        vol1, vol2, vol3 = [], [], []

        vol1_df, vol2_df, vol3_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            if self.__isData1(W, Q2):
                vol1.append([W, Q2])
            if self.__isData2(W, Q2):
                vol2.append([W, Q2])
            if self.__isData3(W, Q2):
                vol3.append([W, Q2])

        if len(vol1) > 0:
            vol1_df = pd.DataFrame(vol1, columns=['W', 'Q2'])
            vol1_df = self.__giveData(vol1_df, 1)
        if len(vol2) > 0:
            vol2_df = pd.DataFrame(vol2, columns=['W', 'Q2'])
            vol2_df = self.__giveData(vol2_df, 2)
        if len(vol3) > 0:
            vol3_df = pd.DataFrame(vol3, columns=['W', 'Q2'])
            vol3_df = self.__giveData(vol3_df, 3)

        result = pd.concat([vol1_df, vol2_df, vol3_df], sort=False)

        result.loc[result['St'] < 0, 'St'] = 0
        result.loc[result['Sl'] < 0, 'Sl'] = 0

        del vol1, vol2, vol3
        del vol1_df, vol2_df, vol3_df

        result = result.groupby(['W', 'Q2', 'cos_th'], as_index=False).mean()

        return result

    def __trig_points_plane(self, Q2_list: list):  # continue
        result = pd.DataFrame()

        buff = self.Data.Data_Q2cos[['cos_th', 'Q2', 'r_TT', 'dr_TT']]

        if self.ratio_str:
            for cos_th in [-1, 1]:
                for Q2 in [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
                    buff = pd.concat([buff, pd.DataFrame({'cos_th': [cos_th], 'Q2': [
                        Q2], 'r_TT': [0], 'dr_TT': [0]})], ignore_index=True)

        x = np.array(buff['cos_th'])
        y = np.array(buff['Q2'])

        f = np.array(buff['r_TT'])
        df = np.array(buff['dr_TT'])

        grid_x, grid_y = np.meshgrid(self.cos_th, Q2_list)

        ratio = interpolate.griddata((x, y), f, (grid_x, grid_y), method='linear')
        dratio = interpolate.griddata((x, y), df, (grid_x, grid_y), method='linear')

        for i in range(grid_x.shape[0]):
            for j in range(grid_y.shape[1]):
                result = pd.concat([result, pd.DataFrame({'Q2': [grid_y[i, j]], 'cos_th': [grid_x[i, j]], 'r_TT': [
                                   ratio[i, j]], 'dr_TT': [dratio[i, j]]})], ignore_index=True)

        del buff, x, y, f, df, grid_x, grid_y, ratio, dratio

        return result

    def __lower_Q2_int(self, data_table):
        result = pd.DataFrame()
        vol1, vol2, vol3, vol4 = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        vol_for_ratio = pd.DataFrame()
        get_ratios = False

        Q2_list = []

        W_Q2_Sigma, W_Q2_None, W_Q2_ratio = {}, {}, {}

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            if self.__isW1(W):
                vol2 = pd.concat(
                    [pd.DataFrame({'W': [W, W], 'Q2':[0.65, 1.0]}), vol2], ignore_index=True)
            elif self.__isW2(W):
                vol2 = pd.concat(
                    [pd.DataFrame({'W': [W, W], 'Q2':[1.0, 1.55]}), vol2], ignore_index=True)
            elif self.__isW3(W):
                vol2 = pd.concat(
                    [pd.DataFrame({'W': [W, W], 'Q2':[1.8, 2.6]}), vol2], ignore_index=True)

            if self.__isSigma(W):
                vol3 = pd.concat([pd.DataFrame({'W': [W], 'Q2':[0.0]}), vol3], ignore_index=True)
                W_Q2_Sigma.setdefault(row['W'], []).append(row['Q2'])
            else:
                if (self.name == "K+L" and self.ratio_str and W > 2.18) or (self.name == "K+S0" and self.ratio_str and W > 2.17):
                    get_ratios = True
                    Q2_list.append(Q2)
                    W_Q2_ratio.setdefault(row['W'], []).append(row['Q2'])
                    vol_for_ratio = pd.concat(
                        [vol_for_ratio, pd.DataFrame({'W': [W], 'Q2': [Q2]})], ignore_index=True)
                else:
                    W_Q2_None.setdefault(row['W'], []).append(row['Q2'])
                    if self.__isW1(W):
                        vol2 = pd.concat(
                            [pd.DataFrame({'W': [W], 'Q2':[0.8]}), vol2], ignore_index=True)
                    elif self.__isW2(W):
                        vol2 = pd.concat(
                            [pd.DataFrame({'W': [W], 'Q2':[1.35]}), vol2], ignore_index=True)
                    elif self.__isW3(W):
                        vol2 = pd.concat(
                            [pd.DataFrame({'W': [W], 'Q2':[2.2]}), vol2], ignore_index=True)
            vol1 = pd.concat([pd.DataFrame({'W': [W], 'Q2':[0.0]}), vol1], ignore_index=True)

        vol1.drop_duplicates(inplace=True)
        vol2.drop_duplicates(inplace=True)
        vol3.drop_duplicates(inplace=True)

        vol1 = self.__Photo_diff(vol1)
        vol2 = self.__Str_func_for_lower(vol2)
        if vol3.shape[0] > 0:
            vol3 = self.__Photo_Sigma(vol3)

            vol3 = vol3.merge(vol1, on=['W', 'Q2', 'cos_th'])
            vol3['Stt'] = -vol3['St']*vol3['Sigma']
            vol3['dStt'] = ((vol3['dSt']*vol3['Sigma'])**2 + (vol3['St']*vol3['dSigma'])**2)**0.5

            vol3.drop(columns=['St', 'dSt', 'Sigma', 'dSigma',
                               'Sl', 'dSl', 'Slt', 'dSlt'], inplace=True)

            vol3 = pd.concat(
                [vol3, vol2.drop(columns=['St', 'dSt', 'Slt', 'dSlt', 'Sl', 'dSl'])], ignore_index=True)

        str3 = pd.DataFrame()
        str3 = pd.concat([vol1, vol2.drop(columns=['Stt', 'dStt'])], ignore_index=True)

        W_Q2_dict = {}

        for index, row in data_table.iterrows():
            W_Q2_dict.setdefault(row['W'], []).append(row['Q2'])

        for W in W_Q2_dict.keys():
            for cos_th in self.cos_th:
                buff = str3.loc[(str3['W'] == W) & (str3['cos_th'] == cos_th)]

                x = np.array(buff['Q2'])

                St = np.array(buff['St'])
                Sl = np.array(buff['Sl'])
                Slt = np.array(buff['Slt'])

                dSt = np.array(buff['dSt'])
                dSl = np.array(buff['dSl'])
                dSlt = np.array(buff['dSlt'])

                f1 = interpolate.interp1d(x, St, kind='quadratic', copy=False)
                f2 = interpolate.interp1d(x, Sl, kind='quadratic', copy=False)
                f3 = interpolate.interp1d(x, Slt, kind='quadratic', copy=False)

                df1 = interpolate.interp1d(x, dSt, kind='linear', copy=False)
                df2 = interpolate.interp1d(x, dSl, kind='linear', copy=False)
                df3 = interpolate.interp1d(x, dSlt, kind='linear', copy=False)

                for Q2 in W_Q2_dict[W]:
                    result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [f1(Q2).astype(np.float)], 'Sl': [f2(Q2).astype(np.float)], 'Slt': [
                        f3(Q2).astype(np.float)], 'dSt': [df1(Q2).astype(np.float)], 'dSl': [df2(Q2).astype(np.float)], 'dSlt': [df3(Q2).astype(np.float)]})], ignore_index=True)

        result_from_Sigma = pd.DataFrame()

        for W in W_Q2_Sigma.keys():
            for cos_th in self.cos_th:
                buff = vol3.loc[(vol3['W'] == W) & (vol3['cos_th'] == cos_th)]

                x = np.array(buff['Q2'])

                Stt = np.array(buff['Stt'])
                dStt = np.array(buff['dStt'])

                f3 = interpolate.interp1d(x, Stt, kind='quadratic', copy=False)
                df3 = interpolate.interp1d(x, dStt, kind='linear', copy=False)

                for Q2 in W_Q2_Sigma[W]:
                    result_from_Sigma = pd.concat([result_from_Sigma, pd.DataFrame(
                        {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'Stt': [f3(Q2).astype(np.float)], 'dStt': [df3(Q2).astype(np.float)]})], ignore_index=True)

        result_from_None = pd.DataFrame()

        for W in W_Q2_None.keys():
            for cos_th in self.cos_th:
                buff = vol2.loc[(vol2['W'] == W) & (vol2['cos_th'] == cos_th)]

                x = np.array(buff['Q2'])

                Stt = np.array(buff['Stt'])
                dStt = np.array(buff['dStt'])

                f3 = interpolate.InterpolatedUnivariateSpline(x, Stt, k=2)
                df3 = interpolate.InterpolatedUnivariateSpline(x, dStt, k=2)

                for Q2 in W_Q2_None[W]:
                    result_from_None = pd.concat([result_from_None, pd.DataFrame(
                        {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'Stt': [f3(Q2).astype(np.float)], 'dStt': [df3(Q2).astype(np.float)]})], ignore_index=True)

        result_from_ratio = pd.DataFrame()

        if vol_for_ratio.shape[0] > 0:
            vol4 = self.__trig_points_plane(Q2_list)
            vol4 = vol_for_ratio.merge(vol4, on='Q2')
            vol4 = result.merge(vol4, on=['W', 'Q2', 'cos_th'])
            vol4['Stt'] = vol4['St']*vol4['r_TT']
            vol4['dStt'] = vol4['dSt']*vol4['r_TT']

            result_from_ratio = (vol4.drop(columns=['St', 'dSt', 'Sl', 'dSl', 'Slt', 'dSlt', 'r_TT', 'dr_TT'])).drop_duplicates(
                subset=['W', 'Q2', 'cos_th'])
            result_from_ratio.reset_index(inplace=True, drop=True)

        result_Stt = pd.DataFrame()

        if(result_from_None.shape[0] > 0):
            result_Stt = pd.concat([result_from_None, result_Stt], ignore_index=True)

        if(result_from_Sigma.shape[0] > 0):
            result_Stt = pd.concat([result_from_Sigma, result_Stt], ignore_index=True)

        if(result_from_ratio.shape[0] > 0):
            result_Stt = pd.concat([result_from_ratio, result_Stt], ignore_index=True)

        result_Stt.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)
        result = result.set_index(['W', 'Q2', 'cos_th'], drop=True)
        result = result.merge(result_Stt, on=['W', 'Q2', 'cos_th'])
        result.reset_index(inplace=True)

        del vol1, vol2, vol3, vol4
        del vol_for_ratio
        del Q2_list, str3, buff
        del result_from_Sigma, result_from_None, result_from_ratio
        del W_Q2_dict, W_Q2_None, W_Q2_Sigma
        del result_Stt

        result['St'] = result['St'].astype(float)
        result['Sl'] = result['Sl'].astype(float)
        result['Slt'] = result['Slt'].astype(float)
        result['Stt'] = result['Stt'].astype(float)
        result['dSt'] = result['dSt'].astype(float)
        result['dSl'] = result['dSl'].astype(float)
        result['dSlt'] = result['dSlt'].astype(float)
        result['dStt'] = result['dStt'].astype(float)

        return result

    def __extrapolate_higher(self, data_table, Q2_list):
        W = data_table.iloc[0]['W']
        cos_th = data_table.iloc[0]['cos_th']

        x = np.array(data_table['Q2'])

        St = data_table['St']
        dSt = data_table['dSt']
        Sl = data_table['Sl']
        dSl = data_table['dSl']
        Slt = data_table['Slt']
        dSlt = data_table['dSlt']
        Stt = data_table['Stt']
        dStt = data_table['dStt']

        def OPE(x, a, b, c):
            y = a + b/x + c/x**2
            return y

        def OPE_p(x, a, b, c):
            y = a**2 + b**2/x + c**2/x**2
            return y

        popt_St, pcov_St = curve_fit(OPE_p, x, St, sigma=dSt, absolute_sigma=True)
        popt_Sl, pcov_Sl = curve_fit(OPE_p, x, Sl, sigma=dSl, absolute_sigma=True)
        popt_Slt, pcov_Slt = curve_fit(OPE, x, Slt, sigma=dSlt, absolute_sigma=True)
        popt_Stt, pcov_Stt = curve_fit(OPE, x, Stt, sigma=dStt, absolute_sigma=True)

        result = pd.DataFrame()

        for Q2 in Q2_list:
            buff = {'W': [W], 'Q2': [Q2], 'cos_th': [cos_th]}

            buff['St'] = [OPE_p(Q2, popt_St[0], popt_St[1], popt_St[2])]
            buff['Sl'] = [OPE_p(Q2, popt_Sl[0], popt_Sl[1], popt_Sl[2])]
            buff['Slt'] = [OPE(Q2, popt_Slt[0], popt_Slt[1], popt_Slt[2])]
            buff['Stt'] = [OPE(Q2, popt_Stt[0], popt_Stt[1], popt_Stt[2])]

            result = pd.concat([result, pd.DataFrame(buff)], ignore_index=True)

        result_for_merge = self.__interp_cub_extraQ2(data_table, Q2_list)
        result_for_merge.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)
        result = result.set_index(['W', 'Q2', 'cos_th'], drop=True)
        result = result.merge(result_for_merge, on=['W', 'Q2', 'cos_th'])
        result.reset_index(inplace=True)

        del x, St, Sl, Stt, Slt, dSt, dSl, dStt, dSlt
        del popt_St, pcov_St, popt_Sl, pcov_Sl
        del popt_Slt, pcov_Slt, popt_Stt, pcov_Stt
        del buff, result_for_merge

        return result

    def __higher_Q2_int(self, data_table):
        result = pd.DataFrame()
        buff = pd.DataFrame()

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            buff = pd.concat([pd.DataFrame([{'W0': W, 'Q20': Q2, 'W': W, 'Q2': 1.8}, {'W0': W, 'Q20': Q2,
                                                                                      'W': W, 'Q2': 2.6}, {'W0': W, 'Q20': Q2, 'W': W, 'Q2': 3.45}]), buff], ignore_index=True)
        W_list = buff['W']
        W_list = W_list.drop_duplicates(keep='first')
        Q2_list = buff['Q20']
        Q2_list = Q2_list.drop_duplicates(keep='first')

        buff = buff.groupby(['W', 'Q2'], as_index=False).mean()
        buff = buff.drop(columns=['W0', 'Q20'])

        buff = self.__giveData(buff, 3)

        for W in W_list:
            for cos_th in self.cos_th:
                result = pd.concat([result, self.__extrapolate_higher(
                    buff.loc[(buff['W'] == W) & (buff['cos_th'] == cos_th)], Q2_list)], ignore_index=True)

        result.loc[result['dSt'] < 0, 'dSt'] = -result.loc[result['dSt'] < 0, 'dSt']
        result.loc[result['dSl'] < 0, 'dSl'] = -result.loc[result['dSl'] < 0, 'dSl']
        result.loc[result['dSlt'] < 0, 'dSlt'] = -result.loc[result['dSlt'] < 0, 'dSlt']
        result.loc[result['dStt'] < 0, 'dStt'] = -result.loc[result['dStt'] < 0, 'dStt']

        del buff, W_list, Q2_list

        return result

    def __Str_func(self, data_table):
        vol1, vol2, vol3, vol4, vol5 = [], [], [], [], []

        vol1_df, vol2_df, vol3_df = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        vol4_df, vol5_df = pd.DataFrame(), pd.DataFrame()

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            if self.__isData1(W, Q2):
                vol1.append([W, Q2])
            if self.__isData2(W, Q2):
                vol2.append([W, Q2])
            if self.__isData3(W, Q2):
                vol3.append([W, Q2])

            if Q2 < 3.45 and self.__isW3(W) and not(self.__isData1(W, Q2) or self.__isData2(W, Q2) or self.__isData3(W, Q2)):
                vol4.append([W, Q2])

            if Q2 > 3.45 and self.__isW3(W) and not(self.__isData1(W, Q2) or self.__isData2(W, Q2) or self.__isData3(W, Q2)):
                vol5.append([W, Q2])

        if len(vol1) > 0:
            vol1_df = pd.DataFrame(vol1, columns=['W', 'Q2'])
            vol1_df = self.__giveData(vol1_df, 1)
        if len(vol2) > 0:
            vol2_df = pd.DataFrame(vol2, columns=['W', 'Q2'])
            vol2_df = self.__giveData(vol2_df, 2)
        if len(vol3) > 0:
            vol3_df = pd.DataFrame(vol3, columns=['W', 'Q2'])
            vol3_df = self.__giveData(vol3_df, 3)
        if len(vol4) > 0:
            vol4_df = pd.DataFrame(vol4, columns=['W', 'Q2'])
            vol4_df = self.__lower_Q2_int(vol4_df)
        if len(vol5) > 0:
            vol5_df = pd.DataFrame(vol5, columns=['W', 'Q2'])
            vol5_df = self.__higher_Q2_int(vol5_df)

        result = pd.concat([vol1_df, vol2_df, vol3_df, vol4_df, vol5_df], ignore_index=True)

        result.loc[result['St'] < 0, 'St'] = 0
        result.loc[result['Sl'] < 0, 'Sl'] = 0

        del vol1, vol2, vol3, vol4, vol5
        del vol1_df, vol2_df, vol3_df, vol4_df, vol5_df

        result = result.groupby(['W', 'Q2', 'cos_th'], as_index=False).mean()

        return result

    def __get_sys_W_err(self, W_list: list):

        result = pd.DataFrame()

        buff = self.Data.W_syserr

        for cos_th in [-1, -0.5, 0, 0.5, 1]:
            buff = pd.concat(
                [buff, pd.DataFrame({'cos_th': [cos_th], 'W': [2.575], 'delta_W': [0]})], ignore_index=True)

        x = np.array(buff['cos_th'])
        y = np.array(buff['W'])

        f = np.array(buff['delta_W'])

        grid_x, grid_y = np.meshgrid(self.cos_th, W_list)

        delta = interpolate.griddata((x, y), f, (grid_x, grid_y), method='linear')

        for i in range(grid_x.shape[0]):
            for j in range(grid_y.shape[1]):
                result = pd.concat([result, pd.DataFrame({'W': [grid_y[i, j]], 'cos_th': [grid_x[i, j]], 'delta_W': [
                                   delta[i, j]]})], ignore_index=True)

        result.loc[result['delta_W'].isnull(), 'delta_W'] = 0.0

        del buff, x, y, f, grid_x, grid_y, delta

        return result

    def __extrapolate_W(self, data_table):
        W_Q2_low, W_Q2_high = {}, {}
        vol1, vol2 = pd.DataFrame(), pd.DataFrame()
        vol3, vol4 = pd.DataFrame(), pd.DataFrame()
        buff = pd.DataFrame()
        result = pd.DataFrame()

        for index, row in data_table.iterrows():
            W = row['W']
            Q2 = row['Q2']

            if W < 2.0:
                W_Q2_low.setdefault(row['Q2'], []).append(row['W'])
                vol1 = pd.concat([vol1, pd.DataFrame({'W': [W], 'Q2': [Q2]})], ignore_index=True)
            else:
                W_Q2_high.setdefault(row['Q2'], []).append(row['W'])
                vol2 = pd.concat([vol2, pd.DataFrame({'W': [W], 'Q2': [Q2]})], ignore_index=True)

        if vol1.shape[0] > 0:
            for Q2 in vol1['Q2'].drop_duplicates():
                if self.name == "K+L":
                    if Q2 < 1.0:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.65], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.725], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.775], 'Q2': [Q2]})], ignore_index=True)
                    elif Q2 >= 1.0 and Q2 < 1.8:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.65], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.75], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.85], 'Q2': [Q2]})], ignore_index=True)
                    elif Q2 >= 1.8:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.63], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.675], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.725], 'Q2': [Q2]})], ignore_index=True)
                else:
                    if Q2 < 1.0:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.725], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.775], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.825], 'Q2': [Q2]})], ignore_index=True)
                    elif Q2 >= 1.0 and Q2 < 1.8:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.75], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.85], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.95], 'Q2': [Q2]})], ignore_index=True)
                    elif Q2 >= 1.8:
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.695], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.725], 'Q2': [Q2]})], ignore_index=True)
                        vol3 = pd.concat(
                            [vol3, pd.DataFrame({'W': [1.775], 'Q2': [Q2]})], ignore_index=True)

        if vol2.shape[0] > 0:
            for Q2 in vol2['Q2'].drop_duplicates():
                vol4 = pd.concat(
                    [vol4, pd.DataFrame({'W': [2.025], 'Q2': [Q2]})], ignore_index=True)
                vol4 = pd.concat(
                    [vol4, pd.DataFrame({'W': [2.275], 'Q2': [Q2]})], ignore_index=True)
                vol4 = pd.concat(
                    [vol4, pd.DataFrame({'W': [2.425], 'Q2': [Q2]})], ignore_index=True)
                vol4 = pd.concat(
                    [vol4, pd.DataFrame({'W': [2.575], 'Q2': [Q2]})], ignore_index=True)

        buff = self.__Str_func(pd.concat([vol3, vol4], ignore_index=True))

        buff.set_index(['W', 'Q2'], inplace=True, drop=True)

        if vol3.shape[0] > 0:
            vol3.set_index(['W', 'Q2'], inplace=True, drop=True)
            vol3 = buff.merge(vol3, on=['W', 'Q2'])
            vol3.reset_index(inplace=True)

        if vol4.shape[0] > 0:
            vol4 = buff.merge(vol4, on=['W', 'Q2'])
            vol4.set_index(['W', 'Q2'], inplace=True, drop=True)
            vol4.reset_index(inplace=True)

        for Q2 in W_Q2_low.keys():
            for cos_th in self.cos_th:
                buff = vol3.loc[(vol3['Q2'] == Q2) & (vol3['cos_th'] == cos_th)]
                if self.name == "K+L":
                    buff = pd.concat([buff, pd.DataFrame({'W': [1.61], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [0], 'dSt': [0],
                                     'Sl': [0], 'dSl': [0], 'Slt': [0], 'dSlt': [0], 'Stt': [0], 'dStt': [0]})], ignore_index=True)
                else:
                    buff = pd.concat([buff, pd.DataFrame({'W': [1.68632], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [0], 'dSt': [0],
                                     'Sl': [0], 'dSl': [0], 'Slt': [0], 'dSlt': [0], 'Stt': [0], 'dStt': [0]})], ignore_index=True)

                x = np.array(buff['W'])

                St = np.array(buff['St'])
                Sl = np.array(buff['Sl'])
                Slt = np.array(buff['Slt'])
                Stt = np.array(buff['Stt'])

                dSt = np.array(buff['dSt'])
                dSl = np.array(buff['dSl'])
                dSlt = np.array(buff['dSlt'])
                dStt = np.array(buff['dStt'])

                f1 = interpolate.interp1d(x, St, kind='quadratic', copy=False)
                f2 = interpolate.interp1d(x, Sl, kind='quadratic', copy=False)
                f3 = interpolate.interp1d(x, Slt, kind='quadratic', copy=False)
                f4 = interpolate.interp1d(x, Stt, kind='quadratic', copy=False)

                df1 = interpolate.interp1d(x, dSt, kind='quadratic', copy=False)
                df2 = interpolate.interp1d(x, dSl, kind='quadratic', copy=False)
                df3 = interpolate.interp1d(x, dSlt, kind='quadratic', copy=False)
                df4 = interpolate.interp1d(x, dStt, kind='quadratic', copy=False)

                for W in W_Q2_low[Q2]:
                    result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [f1(W).astype(np.float)], 'Sl': [f2(W).astype(np.float)], 'Slt': [
                        f3(W).astype(np.float)], 'Stt': [f4(W).astype(np.float)], 'dSt': [df1(W).astype(np.float)], 'dSl': [df2(W).astype(np.float)], 'dSlt': [df3(W).astype(np.float)], 'dStt': [df4(W).astype(np.float)]})], ignore_index=True)

        for Q2 in W_Q2_high.keys():
            for cos_th in self.cos_th:
                buff = vol4.loc[(vol4['Q2'] == Q2) & (vol4['cos_th'] == cos_th)]

                x = np.array(buff['W'])

                St = np.array(buff['St'])
                Sl = np.array(buff['Sl'])
                Slt = np.array(buff['Slt'])
                Stt = np.array(buff['Stt'])

                dSt = np.array(buff['dSt'])
                dSl = np.array(buff['dSl'])
                dSlt = np.array(buff['dSlt'])
                dStt = np.array(buff['dStt'])

                e1 = interpolate.InterpolatedUnivariateSpline(x, St, k=2)
                e2 = interpolate.InterpolatedUnivariateSpline(x, Sl, k=2)
                e3 = interpolate.InterpolatedUnivariateSpline(x, Slt, k=2)
                e4 = interpolate.InterpolatedUnivariateSpline(x, Stt, k=2)

                de1 = interpolate.InterpolatedUnivariateSpline(x, dSt, k=2)
                de2 = interpolate.InterpolatedUnivariateSpline(x, dSl, k=2)
                de3 = interpolate.InterpolatedUnivariateSpline(x, dSlt, k=2)
                de4 = interpolate.InterpolatedUnivariateSpline(x, dStt, k=2)

                for W in W_Q2_high[Q2]:
                    result = pd.concat([result, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'St': [e1(W).astype(np.float)], 'Sl': [e2(W).astype(np.float)], 'Slt': [
                        e3(W).astype(np.float)], 'Stt': [e4(W).astype(np.float)], 'dSt': [de1(W).astype(np.float)], 'dSl': [de2(W).astype(np.float)], 'dSlt': [de3(W).astype(np.float)], 'dStt': [de4(W).astype(np.float)]})], ignore_index=True)

        if self.W_sys and vol2.shape[0] > 0:
            buff = self.__get_sys_W_err(vol2['W'].drop_duplicates())

            buff.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result.set_index(['W', 'cos_th'], inplace=True, drop=True)
            result = pd.merge(buff, result, on=['W', 'cos_th'], how='right')
            result.reset_index(inplace=True)

            result.loc[result['delta_W'].notnull(), 'dSt'] = (
                result['dSt']**2 + result['delta_W']**2)**0.5

            result = result.drop(columns=['delta_W'])

        del buff
        del W_Q2_low, W_Q2_high
        del vol1, vol2
        del vol3, vol4

        return result

    def Str_func_all(self):
        print(' Stand by for structure functions ...')

        vol, vol1, vol2 = [], [], []

        vol_inter, vol_extra, vol_zero = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        for W in self.W:
            if (self.name == "K+L" and W <= 1.61) or (self.name == "K+S0" and W <= 1.68632):
                for Q2 in self.Q2:
                    for cos_th in self.cos_th:
                        vol.append([W, Q2, cos_th])
            else:
                for Q2 in self.Q2:
                    if self.__isWforQ2inter(W, Q2):
                        vol1.append([W, Q2])
                    else:
                        vol2.append([W, Q2])

        if len(vol1) > 0:
            vol_inter = pd.DataFrame(vol1, columns=['W', 'Q2'])
            vol_inter = self.__Str_func(vol_inter)
        if len(vol2) > 0:
            vol_extra = pd.DataFrame(vol2, columns=['W', 'Q2'])
            vol_extra = self.__extrapolate_W(vol_extra)
        if len(vol) > 0:
            vol_zero = pd.DataFrame(vol, columns=['W', 'Q2', 'cos_th'])
            vol_zero['St'] = 0
            vol_zero['dSt'] = 0
            vol_zero['Sl'] = 0
            vol_zero['dSl'] = 0
            vol_zero['Slt'] = 0
            vol_zero['dSlt'] = 0
            vol_zero['Stt'] = 0
            vol_zero['dStt'] = 0

        result = pd.concat([vol_inter, vol_extra, vol_zero], ignore_index=True, sort=True)

        result.loc[result['St'] < 0, 'St'] = 0
        result.loc[result['Sl'] < 0, 'Sl'] = 0

        del vol, vol1, vol2
        del vol_inter, vol_extra, vol_zero

        return result.astype(float)

    def Point_diff(self):
        result = self.req_volume
        result.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)

        S = self.Str_func_all()

        S.set_index(['W', 'Q2', 'cos_th'], inplace=True, drop=True)

        print(' Stand by for cross sections ...')

        result = S.merge(result, on=['W', 'Q2', 'cos_th'])
        result.reset_index(inplace=True)

        result['cs'] = result['St'] + self.eps(result['W'], result['Q2'])*result['Sl'] + self.eps(result['W'], result['Q2'])*result['Stt']*cos(
            2*result['phi']*pi/180) + (self.eps(result['W'], result['Q2'])*(self.eps(result['W'], result['Q2']) + 1))**0.5*result['Slt']*cos(result['phi']*pi/180)
        result['dcs'] = ((result['dSt'])**2 + (self.eps(result['W'], result['Q2'])*result['dSl'])**2 + (self.eps(result['W'], result['Q2'])*result['dStt']*cos(
            2*result['phi']*pi/180))**2 + ((self.eps(result['W'], result['Q2'])*(self.eps(result['W'], result['Q2']) + 1))**0.5*result['dSlt']*cos(result['phi']*pi/180))**2)**0.5

        result.loc[result['cs'] < 0, 'cs'] = 0

        del S

        return result

    @classmethod
    def Average_diff(cls, bin_list: list, options: dict):
        result = pd.DataFrame()
        result_final = pd.DataFrame()

        for bin in bin_list:
            result = pd.concat([result, bin.DataFrame], ignore_index=True)

        config_model = Model(options, result['W'].drop_duplicates().tolist(), result['Q2'].drop_duplicates(
        ).tolist(), result['cos_th'].drop_duplicates().tolist(), result['phi'].drop_duplicates().tolist())
        print(config_model)
        result = config_model.Point_diff()

        for bin in bin_list:
            result_final = pd.concat([result_final, pd.merge(bin.DataFrame, result,
                                                             on=['W', 'Q2', 'cos_th', 'phi']).mean(axis=0).to_frame().T], ignore_index=True)

        del result

        result_final.drop(columns=['St', 'dSt', 'Sl', 'dSl', 'Slt',
                          'dSlt', 'Stt', 'dStt'], inplace=True)
        result_final.rename(columns={'W': 'W_avg', 'Q2': 'Q2_avg', 'cos_th': 'cos_th_avg',
                            'phi': 'phi_avg', 'cs': 'cs_avg', 'dcs': 'dcs_avg'}, inplace=True)

        return result_final

    @classmethod
    def Average_diff_dependency(cls, phase_space: dict, options: dict, W_axis=False, Q2_axis=False, cos_th_axis=False, phi_axis=False):

        bins = []

        if W_axis and Q2_axis:
            for W in range(int(1/phase_space['delta_W'])):
                for Q2 in range(int(5/phase_space['delta_Q2'])):
                    bins.append(Bin({'W': [1.61 + W*phase_space['delta_W'], 1.61 + (W + 1)*phase_space['delta_W']], 'Q2': [Q2*phase_space['delta_Q2'], (Q2 + 1)*phase_space['delta_Q2']], 'cos_th': [phase_space['cos_th'] -
                                phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)

        elif W_axis and cos_th_axis:
            for W in range(int(1/phase_space['delta_W'])):
                for cos_th in range(int(2/phase_space['delta_cos_th'])):
                    bins.append(Bin({'W': [1.61 + W*phase_space['delta_W'], 1.61 + (W + 1)*phase_space['delta_W']], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] + phase_space['delta_Q2']/2],
                                'cos_th': [-1 + cos_th*phase_space['delta_cos_th'], -1 + (cos_th+1)*phase_space['delta_cos_th']], 'phi': [phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)
        elif W_axis and phi_axis:
            for W in range(int(1/phase_space['delta_W'])):
                for phi in range(int(360/phase_space['delta_phi'])):
                    bins.append(Bin({'W': [1.61 + W*phase_space['delta_W'], 1.61 + (W + 1)*phase_space['delta_W']], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] + phase_space['delta_Q2']/2], 'cos_th': [
                        phase_space['cos_th'] - phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phi*phase_space['delta_phi'], (phi+1)*phase_space['delta_phi']]}))
            return Model.Average_diff(bins, options)
        elif Q2_axis and cos_th_axis:
            for Q2 in range(int(5/phase_space['delta_Q2'])):
                for cos_th in range(int(2/phase_space['delta_cos_th'])):
                    bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [Q2*phase_space['delta_Q2'], (Q2 + 1)*phase_space['delta_Q2']], 'cos_th': [-1 +
                                cos_th*phase_space['delta_cos_th'], -1 + (cos_th+1)*phase_space['delta_cos_th']], 'phi': [phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)
        elif Q2_axis and phi_axis:
            for Q2 in range(int(5/phase_space['delta_Q2'])):
                for phi in range(int(360/phase_space['delta_phi'])):
                    bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [Q2*phase_space['delta_Q2'], (Q2 + 1)*phase_space['delta_Q2']], 'cos_th': [
                        phase_space['cos_th'] - phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phi*phase_space['delta_phi'], (phi+1)*phase_space['delta_phi']]}))
            return Model.Average_diff(bins, options)
        elif cos_th_axis and phi_axis:
            for cos_th in range(int(2/phase_space['delta_cos_th'])):
                for phi in range(int(360/phase_space['delta_phi'])):
                    bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] +
                                phase_space['delta_Q2']/2], 'cos_th': [-1 + cos_th*phase_space['delta_cos_th'], -1 + (cos_th+1)*phase_space['delta_cos_th']], 'phi': [phi*phase_space['delta_phi'], (phi+1)*phase_space['delta_phi']]}))
            return Model.Average_diff(bins, options)
        elif W_axis:
            for W in range(int(1/phase_space['delta_W'])):
                bins.append(Bin({'W': [1.61 + W*phase_space['delta_W'], 1.61 + (W + 1)*phase_space['delta_W']], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] + phase_space['delta_Q2']/2], 'cos_th': [
                            phase_space['cos_th'] - phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)
        elif Q2_axis:
            for Q2 in range(int(5/phase_space['delta_Q2'])):
                bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [Q2*phase_space['delta_Q2'], (Q2 + 1)*phase_space['delta_Q2']], 'cos_th': [
                            phase_space['cos_th'] - phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)
        elif cos_th_axis:
            for cos_th in range(int(2/phase_space['delta_cos_th'])):
                bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] + phase_space['delta_Q2']/2], 'cos_th': [-1 + cos_th*phase_space['delta_cos_th'], -1 + (cos_th+1)*phase_space['delta_cos_th']], 'phi': [
                            phase_space['phi'] - phase_space['delta_phi']/2, phase_space['phi'] + phase_space['delta_phi']/2]}))
            return Model.Average_diff(bins, options)
        elif phi_axis:
            for phi in range(int(360/phase_space['delta_phi'])):
                bins.append(Bin({'W': [phase_space['W'] - phase_space['delta_W']/2, phase_space['W'] + phase_space['delta_W']/2], 'Q2': [phase_space['Q2'] - phase_space['delta_Q2']/2, phase_space['Q2'] + phase_space['delta_Q2']/2], 'cos_th': [
                            phase_space['cos_th'] - phase_space['delta_cos_th']/2, phase_space['cos_th'] + phase_space['delta_cos_th']/2], 'phi': [phi*phase_space['delta_phi'], (phi+1)*phase_space['delta_phi']]}))
            return Model.Average_diff(bins, options)

        result = pd.DataFrame()

        return result
