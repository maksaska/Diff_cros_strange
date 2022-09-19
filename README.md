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

Use 'config.yaml' file to set all KY model parameters. You will also need this file to fill your request. There you can find 3 modules:

 **Model**: This module sets the configuration. There are 6 parameters you can change:

  1. _ratio_str_: (recommended True) If True, Stt structure function is estimated from "ratio" method in the "weak" area of phase space. (See more in CLAS12 analysis note(in progress))
  2. _add_factor_: (recommended True) If True, Slt and Stt structure functions approximation fits along cos_th axis are factorized with sin_th and sin^2_th respectivly.
  3. W_sys: (recommended True) If True, adds an additional error to the output in the area of W extrapolation (with W > 2.575 GeV).
  4. err_option: 


  : True
  : True
  : 3
  channel: "K+L"
  E_beam: 6.535
