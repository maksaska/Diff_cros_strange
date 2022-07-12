from model import Model
from bin import Bin
from drawing import plot_cs_W, plot_cs_Q2, plot_cs_cos, plot_cs_phi, plot_str_W, plot_str_cos, plot_str_Q2
from timer import *
import numpy as np

# Timer start
start_time = time.time()

# ------- Single launch --------

#      Model parameters
ratio_str = True           # Метод для Stt из отношений Stt/St
add_factor = True          # Синусы в фитах Slt Stt по cos_th
W_sys = True               # Дополнительная ошибка в экстраполяции St по W (при W > 2.575 GeV)
err_option = 3  # 1, 2, 3  # 1 - константа; 2 - линецная; 3 - квадратичная экстраполяция ошибок в фитах по cos_th

#        Phase space
W = [1.64, 1.61, 1.9, 2.3, 2.6, 2.4]       # Набор W, в которых запрашивается посчитать
Q2 = [0.5, 1.5, 2.0, 0.8, 0.9, 4, 4.5]     # Набор Q2, в которых запрашивается посчитать
cos_th = [-0.95, 0.65, 0.7]                # Набор cos_th, в которых запрашивается посчитать
phi = np.arange(-180, 180, 45)             # Набор phi, в которых запрашивается посчитать
E_beam = 6.535                             # Энергия пучка

#         Model init

config_model = Model("K+L", W, Q2, cos_th, phi, E_beam, ratio_str, add_factor,
                     W_sys, err_option)  # Инициализация программы перед единичным расчетом
# Вывод настроек расчетов
print(config_model)
# Вывод pd.DataFrame() с структурными функциями
print(config_model.Str_func_all())
# Вывод pd.DataFrame() с дифференциальными сечениями
print(config_model.Point_diff())

# --------- Average launch --------

ratio_str = True
add_factor = True
W_sys = True
err_option = 3  # 1, 2, 3

E_beam = 6.535

# Бин-1: Bin([W_min, W_max], кол-во узлов W_n, [Q2_min, Q2_max], Q2_n, [cos_th_min, cos_th_max], cos_th_n, [phi_min, phi_max], phi_n)
# Лист из бинов, где необходимо посчитать среднее дифф сечение
bins = []
bins.append(Bin([1.61, 1.65], 3, [0, 1.5], 3, [-1, 1], 5, [0, 360], 9))             # Бин-1
bins.append(Bin([1.65, 1.69], 3, [0, 1.5], 3, [-1, 1], 5, [0, 360], 9))             # Бин-2
bins.append(Bin([1.69, 1.73], 3, [0, 1.5], 3, [-1, 1], 5, [0, 360], 9))             # Бин-3
bins.append(Bin([1.73, 1.78], 3, [0, 1.5], 3, [-1, 1], 5, [0, 360], 9))             # Бин-4
bins.append(Bin([1.73, 1.78], 10, [1.5, 1.8], 10, [-1, -0.5], 10, [0, 45], 10))     # Бин-5

# average_list - pd.DataFrame() с средними значениями
average_list = Model.Average_diff(bins, E_beam, ratio_str, add_factor, W_sys, err_option)

print(average_list)

# -------------------------------------

# Timer end
end_time = time.time()
time_lapsed = end_time - start_time
time_convert(time_lapsed)


# DRAWING

# print(config_model.Str_func_all())

#plot_cs_W(config_model, 0.8, 0.0, 6.535, 10)
#plot_cs_Q2(config_model, 1.8, 0.0, 6.535, 10)
#plot_cs_cos(config_model, 1.8, 0.5, 6.535, 10)
#plot_cs_phi(config_model, 1.8, 0.5, 0.5, 6.535)

#plot_str_W("St", config_model, 0.8, 0.1, 2.567)
#plot_str_W("Su", config_model, 1.0, 0.35, 2.567)
#plot_str_W("Slt", config_model, 1.0, 0.35, 4.056)
#plot_str_cos("Slt", config_model, 1.725, 0.65, 2.567)
#plot_str_cos("Stt", config_model, 2.375, 0.0, 2.567)

#plot_str_Q2("Stt", config_model, 1.615, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 1.64, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 1.7, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 1.8, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 1.9, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.0, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.1, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.2, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.3, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.4, 0.5, 2.567)
#plot_str_cos("Stt", config_model, 2.375, 0.0, 2.567)
#plot_str_Q2("Stt", config_model, 2.5, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.6, 0.6, 2.567)
#plot_str_Q2("Stt", config_model, 2.65, 0.6, 2.567)

#plot_str_Q2("St", config_model, 2.6, 0.4, 2.567)
#plot_str_Q2("Sl", config_model, 2.6, 0.4, 2.567)
#plot_str_Q2("Slt", config_model, 2.6, 0.4, 2.567)
#plot_str_Q2("Stt", config_model, 2.6, 0.4, 2.567)

#plot_cs_phi(config_model, 1.8, 0.5, 0.5, 6.535)
