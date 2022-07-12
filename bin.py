import pandas as pd
import numpy as np


class Bin:
    def __init__(self, W_list=[1.7, 1.8], N_W=10, Q2_list=[0.5, 1.0], N_Q2=10, cos_th_list=[0, 0.5], N_cos_th=10, phi_list=[0, 45], N_phi=10):
        # Check the boundaries
        assert W_list[0] <= W_list[1], f"W_min ({W_list[0]} should be less than W_max ({W_list[1]})"
        assert Q2_list[0] <= Q2_list[1], f"Q2_min ({Q2_list[0]} should be less than Q2_max ({Q2_list[1]})"
        assert cos_th_list[0] <= cos_th_list[1], f"cos_th_min ({cos_th_list[0]} should be less than cos_th_max ({cos_th_list[1]})"
        assert phi_list[0] <= phi_list[1], f"phi_min ({phi_list[0]} should be less than phi_max ({phi_list[1]})"

        assert cos_th_list[0] <= 1 and cos_th_list[0] >= - \
            1, f"cos_th_min ({cos_th_list[0]} should be within [-1; 1] interval"
        assert cos_th_list[1] <= 1 and cos_th_list[1] >= - \
            1, f"cos_th_max ({cos_th_list[1]} should be within [-1; 1] interval"

        self.__W_min = W_list[0]
        self.__W_max = W_list[1]

        self.__Q2_min = Q2_list[0]
        self.__Q2_max = Q2_list[1]

        self.__cos_th_min = cos_th_list[0]
        self.__cos_th_max = cos_th_list[1]

        self.__phi_min = phi_list[0]
        self.__phi_max = phi_list[1]

        self.__N_W = N_W
        self.__N_Q2 = N_Q2
        self.__N_cos_th = N_cos_th
        self.__N_phi = N_phi

        buff = pd.DataFrame()

        for W in np.linspace(self.__W_min, self.__W_max, num=self.__N_W):
            for Q2 in np.linspace(self.__Q2_min, self.__Q2_max, num=self.__N_Q2):
                for cos_th in np.linspace(self.__cos_th_min, self.__cos_th_max, num=self.__N_cos_th):
                    for phi in np.linspace(self.__phi_min, self.__phi_max, num=self.__N_phi):
                        buff = pd.concat(
                            [buff, pd.DataFrame({'W': [W], 'Q2': [Q2], 'cos_th': [cos_th], 'phi': [phi]})], ignore_index=True)

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
