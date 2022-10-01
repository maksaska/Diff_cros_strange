# DiffCS_strange

This 2.0 version of the newest differential cross section calculator for exclusive single kaon electroproduction allows you to evaluate the differential cross section value for the large invariants' scales and any <a href="https://www.codecogs.com/eqnedit.php?latex=\cos{\theta}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\cos{\theta}" title="\cos{\theta}" /></a> values. Relatively fast and effective procedures make this program a convenient and reliable choice for data analysis in particle physics.

### Formalism 
The differential cross section of the kaons electroproduction off proton in the one-photon approximation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" title="\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega} = \Gamma \dfrac{d\sigma}{d\Omega}_{\gamma^*} \;\;\;\;\;\;\;\; \Gamma = \dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2 - m^2)E_f}{2mE_i}\dfrac{1}{1 - \epsilon}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;=&space;\dfrac{d\sigma_T}{d\Omega}&space;&plus;&space;\varepsilon\dfrac{d\sigma_L}{d\Omega}&space;&plus;&space;\sqrt{\varepsilon(1&plus;\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi}&space;&plus;&space;\varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;=&space;\dfrac{d\sigma_T}{d\Omega}&space;&plus;&space;\varepsilon\dfrac{d\sigma_L}{d\Omega}&space;&plus;&space;\sqrt{\varepsilon(1&plus;\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi}&space;&plus;&space;\varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" title="\dfrac{d\sigma}{d\Omega}_{\gamma^*} = \dfrac{d\sigma_T}{d\Omega} + \varepsilon\dfrac{d\sigma_L}{d\Omega} + \sqrt{\varepsilon(1+\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi} + \varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" /></a>

## Usage

### Installation

Python is used as a primer programming language for that project. To use it you will need the following libraries: numpy, scipy, pandas and plotly. You can check the libriaries with "pip freeze" command. Make sure you have updated all packages (e.g. 'sudo apt update' and 'sudo apt upgrade')

To get this distributive use 'git clone -b python-version https://github.com/Maksaska/Diff_cros_strange.git'

### How to start?

Use 'config.yaml' file to set all KY model parameters. You will also need this file to fill your request. There you can find 4 modules:

**model**: This module sets the configuration. There are 6 parameters you can change:

  1. _ratio_str_: (recommended True) If True, Stt structure function is estimated from "ratio" method in the "weak" area of phase space. (See more in CLAS12 analysis note(in progress))
  2. _add_factor_: (recommended True) If True, Slt and Stt structure functions approximation fits along cos_th axis are factorized with sin_th and sin^2_th respectivly.
  3. _W_sys_: (recommended True) If True, adds an additional error to the output in the area of W extrapolation (with W > 2.575 GeV).
  4. _err_option_: 1 - constant extrapolation of errors along cos_th axis; 2 - linear extrapolation of errors along cos_th axis (up to a 100% of the last known data point error); 3 - quadratic extrapolation of errors along cos_th axis
  5. _channel_: "K+L" or "K+S0"
  6. _E_beam_: beam energy in GeV
  
**cs_grid**: This module provides you with either cross secton or structure function DataFrames in csv format

  1. _str_func_: if True, creates the DataFrame with all structure functions values for requested phase space points in csv format under './Result/str_fun/currrent_date.csv' path. Phase space grid is formed from all the possible combinations of (W, Q2, cos_th) provided by the user.
  2. _diff_cs_: if True, creates the DataFrame with all cross section values for requested phase space points in csv format under './Result/diff_cs/currrent_date.csv' path. Phase space grid is formed from all the possible combinations of (W, Q2, cos_th, phi) provided by the user.
  3. _W_: [1.7, 1.8]
  4. _Q2_: [0.5, 0.8]
  5. _cos_th_: [0.2, 0.6]
  6. _phi_: [-20, 30] 

**av_cs_grid**: This module provides you with average cross secton DataFrames in csv format. There are 2 methods implemented:

  1. Manual bin input. The user can enter any amount of bins of any size. To add bin, add the following line into module:
    
    binN: {'W': [W_min, W_max], 'Q2': [Q2_min, Q2_max], 'cos_th': [cos_th_min, cos_th_max], 'phi': [phi_min, phi_max]}
    
* To activate the calculations set 'active' input field as 'True' 
    
* The output will be stored in csv format under './Result/av_diff_cs/currrent_date.csv' path. Phase space coordinats in the csv files are average points of specified bins.
    
 2. Average cross section calculation for predetermined number of bin along one (or two) of the axes. Set average phase space point with (W, Q2, cos_th, phi) and bin sizes with (dW, dQ2, dcos_th, dphi). By choosing one (or two) of the axes, you specify which reaction parameter will vary. 
  
  For example: 
      
      
      
      method2:
         active: True
         plot: False
         W_axis: False
         Q2_axis: True
         cos_th_axis: True
         phi_axis: False
         W: 2.0
         Q2: 1.5
         cos_th: 0
         phi: 45
         dW: 0.1
         dQ2: 0.5
         dcos_th: 0.5
         dphi: 10

* _method2_ configuration above will request a DataFrame with average cross section dependencies along Q2 and cos_th with fixed average W and phi variables. Bin sizes for W and phi axes are fixed and equal to dW and dphi. This program will do the calculation for all possible phase space (e.g. for Q2 from 0 to 5 GeV2; W from MY to 2.65 GeV; cos_th from -1 to 1; phi from 0 to 360) Number of points for Q2 and cos_th axes depends on dQ2 and dcos_th. In our case, with 'dcos_th: 0.5' we will get 5 points along cos_th axis: -1, -0.5, 0, 0.5, 1.0. The same goes for Q2 with 'dQ2: 0.5', which runs from 0 to 5 GeV2.
    
  
* To activate the calculations set 'active' input field as 'True' 

* if 'plot' is set to True, provide the user with interactive plot in html format under'./Result/Plots/observable_currrent_date.html' path.

* The output will be stored in csv format under './Result/av_diff_cs/currrent_date.csv' path. Phase space coordinats in the csv files are average points of specified bins.

**plot**: This module provides you with plots for cross section or structure function.

  1. _plot_cs_ - if set to True, creates a plot for cross section for specified (W, Q2, cos_th, phi)
  2. _plot_str_func_ - if set to True, creates a plot for structure functions for specified (W, Q2, cos_th)
  
* To activate the calculations set 'active' input field as 'True' 

* The output will be stored in html format under './Result/Plots/observable_currrent_date.html' path. **Note**: you can only choose one of the axes.
  

  
