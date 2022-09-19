from model import Model
from bin import Bin
from timer import *
from drawing import *
from datetime import datetime
import numpy as np
import plotly.express as px
import yaml
import os


def main():
    with open('config.yaml', 'r') as file:
        configuration = yaml.safe_load(file)

    start_time = time.time()
    today = datetime.now()
    tag = today.strftime("%b_%d_%Y_%H_%M_%S")

    os.makedirs('./Results/str_fun', exist_ok=True)
    os.makedirs('./Results/diff_cs', exist_ok=True)
    os.makedirs('./Results/av_diff_cs', exist_ok=True)
    os.makedirs('./Results/Plots', exist_ok=True)

    options = {'ratio_str': configuration['model']['ratio_str'], 'add_factor': configuration['model']['add_factor'], 'W_sys': configuration['model']
               ['W_sys'], 'err_option': configuration['model']['err_option'], 'E_beam': configuration['model']['E_beam'], 'channel': configuration['model']['channel']}

    # ------- Single launch --------

    if configuration['cs_grid']['str_func'] + configuration['cs_grid']['diff_cs']:
        W = configuration['cs_grid']['W']
        Q2 = configuration['cs_grid']['Q2']
        cos_th = configuration['cs_grid']['cos_th']
        phi = configuration['cs_grid']['phi']

        config_model = Model(options, W, Q2, cos_th, phi)

        if configuration['cs_grid']['str_func']:
            config_model.Str_func_all().to_csv(f'./Results/str_fun/{tag}.csv', index=False)

        if configuration['cs_grid']['diff_cs']:
            config_model.Point_diff().to_csv(f'./Results/diff_cs/{tag}.csv', index=False)

    # --------- Average cross section (method_1) --------

    if configuration['av_cs_grid']['method1']['active']:
        bins = []
        for i in configuration['av_cs_grid']['method1']:
            if i == 'active':
                continue
            bins.append(Bin(configuration['av_cs_grid']['method1'][i]))

        Model.Average_diff(bins, options).to_csv(f'./Results/av_diff_cs/{tag}.csv', index=False)

    # --------- Average cross section (method_2) --------

    if configuration['av_cs_grid']['method2']['active']:
        phase_space = {'W': configuration['av_cs_grid']['method2']['W'], 'Q2': configuration['av_cs_grid']['method2']['Q2'], 'cos_th': configuration['av_cs_grid']['method2']['cos_th'], 'phi': configuration['av_cs_grid']['method2']['phi'],
                       'delta_W': configuration['av_cs_grid']['method2']['dW'], 'delta_Q2': configuration['av_cs_grid']['method2']['dQ2'], 'delta_cos_th': configuration['av_cs_grid']['method2']['dcos_th'], 'delta_phi': configuration['av_cs_grid']['method2']['dphi']}

        Model.Average_diff_dependency(phase_space, options, W_axis=configuration['av_cs_grid']['method2']['W_axis'], Q2_axis=configuration['av_cs_grid'][
                                      'method2']['Q2_axis'], cos_th_axis=configuration['av_cs_grid']['method2']['cos_th_axis'], phi_axis=configuration['av_cs_grid']['method2']['phi_axis']).to_csv(f'./Results/av_diff_cs/{tag}.csv', index=False)

        if configuration['av_cs_grid']['method2']['plot']:
            Show_average_cs(phase_space, options,  W_axis=configuration['av_cs_grid']['method2']['W_axis'], Q2_axis=configuration['av_cs_grid']['method2']
                            ['Q2_axis'], cos_th_axis=configuration['av_cs_grid']['method2']['cos_th_axis'], phi_axis=configuration['av_cs_grid']['method2']['phi_axis']).write_html(f'./Results/Plots/av_cs_{tag}.html', include_mathjax='cdn')

    # --------- Plots --------

    if configuration['plot']['active']:
        phase_space = {'W': [configuration['plot']['W']], 'Q2': [configuration['plot']['Q2']], 'cos_th': [
            configuration['plot']['cos_th']], 'phi': [configuration['plot']['phi']]}

        if configuration['plot']['plot_cs']:
            Show_cs(phase_space, options, W_axis=configuration['plot']['W_axis'], Q2_axis=configuration['plot']
                    ['Q2_axis'], cos_th_axis=configuration['plot']['cos_th_axis'], phi_axis=configuration['plot']['phi_axis']).write_html(f'./Results/Plots/cs_{tag}.html', include_mathjax='cdn')

        if configuration['plot']['plot_str_func']:
            Show_str_func(phase_space, options, W_axis=configuration['plot']['W_axis'], Q2_axis=configuration['plot']
                          ['Q2_axis'], cos_th_axis=configuration['plot']['cos_th_axis']).write_html(f'./Results/Plots/str_func_{tag}.html', include_mathjax='cdn')

    # -------------------------------------

    # Timer end
    end_time = time.time()
    time_lapsed = end_time - start_time
    time_convert(time_lapsed)


if __name__ == "__main__":
    main()
