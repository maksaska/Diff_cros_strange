import pandas as pd
import numpy as np


class Bin:
    def __init__(self, ps_dict):
        # Check the boundaries
        assert ps_dict['W'][0] <= ps_dict['W'][1], f"W_min ({ps_dict['W'][0]} should be less than W_max ({ps_dict['W'][1]})"
        assert ps_dict['Q2'][0] <= ps_dict['Q2'][1], f"Q2_min ({ps_dict['Q2'][0]} should be less than Q2_max ({ps_dict['Q2'][1]})"
        assert ps_dict['cos_th'][0] <= ps_dict['cos_th'][
            1], f"cos_th_min ({ps_dict['cos_th'][0]} should be less than cos_th_max ({ps_dict['cos_th'][1]})"
        assert ps_dict['phi'][0] <= ps_dict['phi'][
            1], f"phi_min ({ps_dict['phi'][0]} should be less than phi_max ({ps_dict['phi'][1]})"

        if ps_dict['cos_th'][0] > 1:
            ps_dict['cos_th'][0] = 1

        if ps_dict['cos_th'][1] > 1:
            ps_dict['cos_th'][1] = 1

        if ps_dict['cos_th'][0] < -1:
            ps_dict['cos_th'][0] = -1

        if ps_dict['cos_th'][1] < -1:
            ps_dict['cos_th'][1] = -1

        self.__W_min = ps_dict['W'][0]
        self.__W_max = ps_dict['W'][1]

        self.__Q2_min = ps_dict['Q2'][0]
        self.__Q2_max = ps_dict['Q2'][1]

        self.__cos_th_min = ps_dict['cos_th'][0]
        self.__cos_th_max = ps_dict['cos_th'][1]

        self.__phi_min = ps_dict['phi'][0]
        self.__phi_max = ps_dict['phi'][1]

        self.__N_W = round((self.__W_max - self.__W_min)/0.01) + 1
        self.__N_Q2 = round((self.__Q2_max - self.__Q2_min)/0.1) + 1
        self.__N_cos_th = round((self.__cos_th_max - self.__cos_th_min)/0.1) + 1
        self.__N_phi = round((self.__phi_max - self.__phi_min)/15) + 1

        buff = pd.DataFrame()

        for W in np.linspace(self.__W_min, self.__W_max, num=self.__N_W):
            for Q2 in np.linspace(self.__Q2_min, self.__Q2_max, num=self.__N_Q2):
                for cos_th in np.linspace(self.__cos_th_min, self.__cos_th_max, num=self.__N_cos_th):
                    for phi in np.linspace(self.__phi_min, self.__phi_max, num=self.__N_phi):
                        buff = pd.concat(
                            [buff, pd.DataFrame({'W': [round(W, 5)], 'Q2': [round(Q2, 5)], 'cos_th': [round(cos_th, 5)], 'phi': [round(phi, 5)]})], ignore_index=True)

        self.__DataFrame = buff

        del buff

    @ property
    def DataFrame(self):
        return self.__DataFrame

    def __repr__(self):
        return f'Bin(["{self.__W_min}", "{self.__W_max}"], "{self.__N_W}", ["{self.__Q2_min}", "{self.__Q2_max}"], "{self.__N_Q2}", ["{self.__cos_th_min}", "{self.__cos_th_max}"], "{self.__N_cos_th}", ["{self.__phi_min}", "{self.__phi_max}"], "{self.__N_phi}")'

    def __str__(self):
        print(self.DataFrame)
        return ' '
