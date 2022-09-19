from model import Model
from bin import Bin
import plotly.express as px
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def Show_average_cs(phase_space: dict, options: dict, W_axis=False, Q2_axis=False, cos_th_axis=False, phi_axis=False):
    if W_axis and Q2_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, W_axis=True, Q2_axis=True)

        title = "cos_th = " + str(phase_space['cos_th']) + " +- " + str(phase_space['delta_cos_th']) + " ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter_3d(average_list, x='W_avg', y='Q2_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='W [GeV]', yaxis_title='Q2 [GeV2]',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif W_axis and cos_th_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, W_axis=True, cos_th_axis=True)

        title = "Q2 = " + str(phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + " GeV2 ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter_3d(average_list, x='W_avg', y='cos_th_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='W [GeV]', yaxis_title='cos_th',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif W_axis and phi_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, W_axis=True, phi_axis=True)

        title = "Q2 = " + str(phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + \
            " GeV2 ; cos_th = " + str(phase_space['cos_th']) + \
            " +- " + str(phase_space['delta_cos_th'])
        fig = px.scatter_3d(average_list, x='W_avg', y='phi_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='W [GeV]', yaxis_title='phi [degree]',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif Q2_axis and cos_th_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, Q2_axis=True, cos_th_axis=True)

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter_3d(average_list, x='Q2_avg', y='cos_th_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='Q2 [GeV2]', yaxis_title='cos_th',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif Q2_axis and phi_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, Q2_axis=True, phi_axis=True)

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; cos_th = " + str(
            phase_space['cos_th']) + " +- " + str(phase_space['delta_cos_th'])
        fig = px.scatter_3d(average_list, x='Q2_avg', y='phi_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='Q2 [GeV2]', yaxis_title='phi [degree]',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif cos_th_axis and phi_axis:
        average_list = Model.Average_diff_dependency(
            phase_space, options, cos_th_axis=True, phi_axis=True)

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; Q2 = " + str(
            phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + " GeV2"
        fig = px.scatter_3d(average_list, x='cos_th_avg', y='phi_avg',
                            z='cs_avg', color='dcs_avg', title=title)
        fig.update_layout(scene=dict(xaxis_title='cos_th', yaxis_title='phi [degree]',
                                     zaxis_title='Diff. cross section [nb]'), coloraxis_colorbar=dict(
            title="Absolute error [nb]"))

    elif W_axis:
        average_list = Model.Average_diff_dependency(phase_space, options, W_axis=True)
        average_list['dW'] = phase_space['delta_W']/2

        title = "Q2 = " + str(phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + " GeV2 ; cos_th = " + str(phase_space['cos_th']) + " +- " + str(phase_space['delta_cos_th']) + " ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter(average_list, x='W_avg', y='cs_avg', error_x='dW', error_y='dcs_avg', labels=dict(
            W_avg='W [GeV]', cs_avg='Diff. cross section [nb]'), title=title)

    elif Q2_axis:
        average_list = Model.Average_diff_dependency(phase_space, options, Q2_axis=True)
        average_list['dQ2'] = phase_space['delta_Q2']/2

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; cos_th = " + str(phase_space['cos_th']) + " +- " + str(phase_space['delta_cos_th']) + " ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter(average_list, x='Q2_avg', y='cs_avg',
                         error_x='dQ2', error_y='dcs_avg', labels=dict(Q2_avg='Q2 [GeV2]', cs_avg='Diff. cross section [nb]'), title=title)

    elif cos_th_axis:
        average_list = Model.Average_diff_dependency(phase_space, options, cos_th_axis=True)
        average_list['dcos_th'] = phase_space['delta_cos_th']/2

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; Q2 = " + str(phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + " GeV2 ; phi = " + str(
            phase_space['phi']) + " +- " + str(phase_space['delta_phi']) + " degree"
        fig = px.scatter(average_list, x='cos_th_avg', y='cs_avg',
                         error_x='dcos_th', error_y='dcs_avg', labels=dict(cos_th_avg='cos_th', cs_avg='Diff. cross section [nb]'), title=title)

    elif phi_axis:
        average_list = Model.Average_diff_dependency(phase_space, options, phi_axis=True)
        average_list['dphi'] = phase_space['delta_phi']/2

        title = "W = " + str(phase_space['W']) + " +- " + str(phase_space['delta_W']) + " GeV ; Q2 = " + str(phase_space['Q2']) + " +- " + str(phase_space['delta_Q2']) + " GeV2 ; cos_th = " + str(
            phase_space['cos_th']) + " +- " + str(phase_space['delta_cos_th'])
        fig = px.scatter(average_list, x='phi_avg', y='cs_avg',
                         error_x='dphi', error_y='dcs_avg', labels=dict(phi_avg='phi [degree]', cs_avg='Diff. cross section [nb]'), title=title)

    # fig.show()
    return fig


def Show_cs(phase_space: dict, options: dict, W_axis=False, Q2_axis=False, cos_th_axis=False, phi_axis=False):
    if W_axis:
        W = np.arange(1.61, 2.65, 0.01)
        cs_model = Model(options, W, phase_space['Q2'], phase_space['cos_th'], phase_space['phi'])

        cs_list = cs_model.Point_diff().sort_values(by=['W']).reset_index(drop=True)

        fig = make_subplots(rows=1, cols=1)

        W = cs_list['W']
        cs = cs_list['cs']
        dcs = cs_list['dcs']

        fig.add_trace(go.Scatter(x=W, y=cs, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=W, y=cs+dcs, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=W, y=cs-dcs, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.update_xaxes(title_text=r'$\text{W [GeV]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\gamma^{*}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)

        title = r'$\text{Diff. cross section } Q^{2} = %.2f \text{ GeV}^{2}; \cos{\theta} = %.2f; \phi = %.1f \text{ degree}$' % (
            phase_space['Q2'][0], phase_space['cos_th'][0], phase_space['phi'][0])

        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    elif Q2_axis:
        Q2 = np.arange(0, 5.0, 0.05)
        cs_model = Model(options, phase_space['W'], Q2, phase_space['cos_th'], phase_space['phi'])

        cs_list = cs_model.Point_diff().sort_values(by=['Q2']).reset_index(drop=True)

        fig = make_subplots(rows=1, cols=1)

        Q2 = cs_list['Q2']
        cs = cs_list['cs']
        dcs = cs_list['dcs']

        fig.add_trace(go.Scatter(x=Q2, y=cs, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=Q2, y=cs+dcs, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=Q2, y=cs-dcs, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.update_xaxes(title_text=r'$Q^{2} \text{ [GeV}^{2}\text{]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\gamma^{*}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)

        title = r'$\text{Diff. cross section } W = %.2f \text{ GeV}; \cos{\theta} = %.2f; \phi = %.1f \text{ degree}$' % (
            phase_space['W'][0], phase_space['cos_th'][0], phase_space['phi'][0])

        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    elif cos_th_axis:
        cos_th = np.arange(-1.0, 1.0, 0.01)
        cs_model = Model(options, phase_space['W'], phase_space['Q2'], cos_th, phase_space['phi'])

        cs_list = cs_model.Point_diff().sort_values(by=['cos_th']).reset_index(drop=True)

        fig = make_subplots(rows=1, cols=1)

        cos_th = cs_list['cos_th']
        cs = cs_list['cs']
        dcs = cs_list['dcs']

        fig.add_trace(go.Scatter(x=cos_th, y=cs, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=cos_th, y=cs+dcs, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=cos_th, y=cs-dcs, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.update_xaxes(title_text=r'$\cos{\theta}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\gamma^{*}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)

        title = r'$\text{Diff. cross section } W = %.2f \text{ GeV}; Q^{2} = %.2f \text{ GeV}^2; \phi = %.1f \text{ degree}$' % (
            phase_space['W'][0], phase_space['Q2'][0], phase_space['phi'][0])

        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    elif phi_axis:
        phi = np.arange(-180, 180, 1)
        cs_model = Model(options, phase_space['W'], phase_space['Q2'], phase_space['cos_th'], phi)

        cs_list = cs_model.Point_diff().sort_values(by=['phi']).reset_index(drop=True)

        fig = make_subplots(rows=1, cols=1)

        phi = cs_list['phi']
        cs = cs_list['cs']
        dcs = cs_list['dcs']

        fig.add_trace(go.Scatter(x=phi, y=cs, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=phi, y=cs+dcs, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=phi, y=cs-dcs, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.update_xaxes(title_text=r'$\phi \text{ [degree]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\gamma^{*}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)

        title = r'$\text{Diff. cross section } W = %.2f \text{ GeV}; Q^{2} = %.2f \text{ GeV}^2; \cos{\theta} = %.2f$' % (
            phase_space['W'][0], phase_space['Q2'][0], phase_space['phi'][0])

        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    # fig.show()
    return fig


def Show_str_func(phase_space: dict, options: dict, W_axis=False, Q2_axis=False, cos_th_axis=False):
    if W_axis:
        W = np.arange(1.61, 2.65, 0.01)
        cs_model = Model(options, W, phase_space['Q2'], phase_space['cos_th'], phase_space['phi'])
        cs_list = cs_model.Str_func_all().sort_values(by=['W']).reset_index(drop=True)

        fig = make_subplots(rows=2, cols=2)

        W = cs_list['W']

        St = cs_list['St']
        Sl = cs_list['Sl']
        Slt = cs_list['Slt']
        Stt = cs_list['Stt']

        dSt = cs_list['dSt']
        dSl = cs_list['dSl']
        dSlt = cs_list['dSlt']
        dStt = cs_list['dStt']

        fig.add_trace(go.Scatter(x=W, y=St, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=W, y=St+dSt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=W, y=St-dSt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.add_trace(go.Scatter(x=W, y=Sl, name="Longitudinal", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=W, y=Sl+dSl, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=W, y=Sl-dSl, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=2)

        fig.add_trace(go.Scatter(x=W, y=Slt, name="LT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=W, y=Slt+dSlt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=W, y=Slt-dSlt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=1)

        fig.add_trace(go.Scatter(x=W, y=Stt, name="TT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=W, y=Stt+dStt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=W, y=Stt-dStt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=2)

        fig.update_xaxes(title_text=r'$\text{W [GeV]}$', row=1, col=1)
        fig.update_xaxes(title_text=r'$\text{W [GeV]}$', row=1, col=2)
        fig.update_xaxes(title_text=r'$\text{W [GeV]}$', row=2, col=1)
        fig.update_xaxes(title_text=r'$\text{W [GeV]}$', row=2, col=2)

        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{T}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{L}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=2)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{LT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{TT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=2)

        title = r'$\text{Structure functions for } Q^{2} = %.2f \text{ GeV }^{2}; \cos{\theta} = %.2f$' % (
            cs_list['Q2'][0], cs_list['cos_th'][0])
        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    elif Q2_axis:
        Q2 = np.arange(0, 5, 0.05)
        cs_model = Model(options, phase_space['W'], Q2, phase_space['cos_th'], phase_space['phi'])
        cs_list = cs_model.Str_func_all().sort_values(by=['Q2'])

        fig = make_subplots(rows=2, cols=2)

        Q2 = cs_list['Q2']

        St = cs_list['St']
        Sl = cs_list['Sl']
        Slt = cs_list['Slt']
        Stt = cs_list['Stt']

        dSt = cs_list['dSt']
        dSl = cs_list['dSl']
        dSlt = cs_list['dSlt']
        dStt = cs_list['dStt']

        fig.add_trace(go.Scatter(x=Q2, y=St, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=Q2, y=St+dSt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=Q2, y=St-dSt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.add_trace(go.Scatter(x=Q2, y=Sl, name="Longitudinal", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=Q2, y=Sl+dSl, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=Q2, y=Sl-dSl, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=2)

        fig.add_trace(go.Scatter(x=Q2, y=Slt, name="LT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=Q2, y=Slt+dSlt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=Q2, y=Slt-dSlt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=1)

        fig.add_trace(go.Scatter(x=Q2, y=Stt, name="TT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=Q2, y=Stt+dStt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=Q2, y=Stt-dStt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=2)

        fig.update_xaxes(title_text=r'$Q^{2} \text{ [GeV}^{2}\text{]}$', row=1, col=1)
        fig.update_xaxes(title_text=r'$Q^{2} \text{ [GeV}^{2}\text{]}$', row=1, col=2)
        fig.update_xaxes(title_text=r'$Q^{2} \text{ [GeV}^{2}\text{]}$', row=2, col=1)
        fig.update_xaxes(title_text=r'$Q^{2} \text{ [GeV}^{2}\text{]}$', row=2, col=2)

        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{T}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{L}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=2)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{LT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{TT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=2)

        title = r'$\text{Structure functions for } W = %.2f \text{ GeV}; \cos{\theta} = %.2f$' % (
            cs_list['W'][0], cs_list['cos_th'][0])
        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    elif cos_th_axis:
        cos_th = np.arange(-1, 1, 0.01)
        cs_model = Model(options, phase_space['W'], phase_space['Q2'], cos_th, phase_space['phi'])
        cs_list = cs_model.Str_func_all().sort_values(by=['cos_th'])

        fig = make_subplots(rows=2, cols=2)

        cos_th = cs_list['cos_th']

        St = cs_list['St']
        Sl = cs_list['Sl']
        Slt = cs_list['Slt']
        Stt = cs_list['Stt']

        dSt = cs_list['dSt']
        dSl = cs_list['dSl']
        dSlt = cs_list['dSlt']
        dStt = cs_list['dStt']

        fig.add_trace(go.Scatter(x=cos_th, y=St, name="Transverse", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=cos_th, y=St+dSt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=cos_th, y=St-dSt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=1)

        fig.add_trace(go.Scatter(x=cos_th, y=Sl, name="Longitudinal", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=1, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=cos_th, y=Sl+dSl, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=cos_th, y=Sl-dSl, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=1, col=2)

        fig.add_trace(go.Scatter(x=cos_th, y=Slt, name="LT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=1)
        fig.add_trace(go.Scatter(name='Upper Bound', x=cos_th, y=Slt+dSlt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=1)
        fig.add_trace(go.Scatter(name='Lower Bound', x=cos_th, y=Slt-dSlt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=1)

        fig.add_trace(go.Scatter(x=cos_th, y=Stt, name="TT", mode='lines',
                      line=dict(color='rgb(31, 119, 180)')), row=2, col=2)
        fig.add_trace(go.Scatter(name='Upper Bound', x=cos_th, y=Stt+dStt, mode='lines',
                      marker=dict(color="#444"), line=dict(width=0), showlegend=False), row=2, col=2)
        fig.add_trace(go.Scatter(name='Lower Bound', x=cos_th, y=Stt-dStt, mode='lines', marker=dict(color="#444"),
                      line=dict(width=0), fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False), row=2, col=2)

        fig.update_xaxes(title_text=r'$\cos{\theta}$', row=1, col=1)
        fig.update_xaxes(title_text=r'$\cos{\theta}$', row=1, col=2)
        fig.update_xaxes(title_text=r'$\cos{\theta}$', row=2, col=1)
        fig.update_xaxes(title_text=r'$\cos{\theta}$', row=2, col=2)

        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{T}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{L}}}{\text{d}\Omega} \text{ [nb]}$', row=1, col=2)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{LT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=1)
        fig.update_yaxes(
            title_text=r'$\dfrac{\text{d}\sigma_{\text{TT}}}{\text{d}\Omega} \text{ [nb]}$', row=2, col=2)

        title = r'$\text{Structure functions for } W = %.2f \text{ GeV}; Q^{2} = %.2f \text{ GeV}^{2}$' % (
            cs_list['W'][0], cs_list['Q2'][0])
        fig.update_layout(height=768, width=1366, hovermode="x", title_text=title)

    # fig.show()
    return fig
