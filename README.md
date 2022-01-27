# DiffCS_strange

This 3.0.1 version of the newest differential cross section calculator for exclusive single kaon electroproduction allows you to evaluate the differential cross section value for the large invariants' scales and any <a href="https://www.codecogs.com/eqnedit.php?latex=\cos{\theta}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\cos{\theta}" title="\cos{\theta}" /></a> values. Relatively fast and effective procedures make this program a convenient and reliable choice for data analysis in particle physics.

### Formalism 
The differential cross section of the kaons electroproduction off proton in the one-photon approximation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" title="\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega} = \Gamma \dfrac{d\sigma}{d\Omega}_{\gamma^*} \;\;\;\;\;\;\;\; \Gamma = \dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2 - m^2)E_f}{2mE_i}\dfrac{1}{1 - \epsilon}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;=&space;\dfrac{d\sigma_T}{d\Omega}&space;&plus;&space;\varepsilon\dfrac{d\sigma_L}{d\Omega}&space;&plus;&space;\sqrt{\varepsilon(1&plus;\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi}&space;&plus;&space;\varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{d\Omega}_{\gamma^*}&space;=&space;\dfrac{d\sigma_T}{d\Omega}&space;&plus;&space;\varepsilon\dfrac{d\sigma_L}{d\Omega}&space;&plus;&space;\sqrt{\varepsilon(1&plus;\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi}&space;&plus;&space;\varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" title="\dfrac{d\sigma}{d\Omega}_{\gamma^*} = \dfrac{d\sigma_T}{d\Omega} + \varepsilon\dfrac{d\sigma_L}{d\Omega} + \sqrt{\varepsilon(1+\varepsilon)}\dfrac{d\sigma_{LT}}{d\Omega}\cos{\varphi} + \varepsilon\dfrac{d\sigma_{TT}}{d\Omega}\cos{2\varphi}" /></a>

## Usage
1. Install [Root Cern](https://root.cern.ch/building-root)
2. Git the program: git clone https://github.com/Maksaska/Diff_cros_strange.git
3. Type command: chmod +x Run
4. Compile with "Run",i.e. ./Run
5. Start the compiled file with ./start

Requirements: [Root Cern](https://root.cern/)

## Options for program start:
* The energy of the incident electron:
You should add "-e X", where X in GeV, to the "./start" command to specify the beam energy.
* The program is written for two channels. The default is <a href="https://www.codecogs.com/eqnedit.php?latex=K\Lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K\Lambda" title="K\Lambda" /></a>. That is if you don't take any additional action.
To work with <a href="https://www.codecogs.com/eqnedit.php?latex=K\Sigma^0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K\Sigma^0" title="K\Sigma^0" /></a> channel, you should add "--KSigma0" to "./start" command like "./start --KSigma0"

* 2 modes of operation are initialized in the program:
  * Cross section calculation at the given point <a href="https://www.codecogs.com/eqnedit.php?latex=(W,&space;Q^2,&space;\cos{\theta},&space;\varphi)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(W,&space;Q^2,&space;\cos{\theta},&space;\varphi)" title="(W, Q^2, \cos{\theta}, \varphi)" /></a>
  * Average cross section calculation for a given area of <a href="https://www.codecogs.com/eqnedit.php?latex=W" target="_blank"><img src="https://latex.codecogs.com/gif.latex?W" title="W" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=Q^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q^2" title="Q^2" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\cos{\theta}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\cos{\theta}" title="\cos{\theta}" /></a>, and <a href="https://www.codecogs.com/eqnedit.php?latex=\varphi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varphi" title="\varphi" /></a>. 

## Average calculation
For this mode of operation, the program requires the boundaries values of the <a href="https://www.codecogs.com/eqnedit.php?latex=W,\;Q^2,\;\cos{\theta},\;\varphi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?W,\;Q^2,\;\cos{\theta},\;\varphi" title="W,\;Q^2,\;\cos{\theta},\;\varphi" /></a> intervals. They are all zero by default.
* --W_min=X1, where X1 is an invariant mass of the final hadron system in <a href="https://www.codecogs.com/eqnedit.php?latex=GeV" target="_blank"><img src="https://latex.codecogs.com/gif.latex?GeV" title="GeV" /></a>
* --W_max=X2
* --Q2_min=X3, where X3 is photon virtuality in <a href="https://www.codecogs.com/eqnedit.php?latex=GeV^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?GeV^2" title="GeV^2" /></a>
* --Q2_max=X4
* --cos_min=X5,where X5 is the cosine value of the kaon polar angle in c.m. frame
* --cos_max=X6
* --phi_min=X7, where X7 is the value of the kaon azimuth angle in c.m. frame
* --phi_max=X8

> An example: ./start -e 6.5 --W_min=1.75 --W_max=1.95 --Q2_min=0.65 --Q2_max=0.7 --cos_min=-1 --cos_max=1 --phi_min=15 --phi_max=45 --KSigma0
> 
> Meaning: <a href="https://www.codecogs.com/eqnedit.php?latex=K\Sigma^0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K\Sigma^0" title="K\Sigma^0" /></a> channel, average cross section calculation, 
> 
> <a href="https://www.codecogs.com/eqnedit.php?latex=E&space;=&space;6.5&space;\text{&space;GeV},&space;\;W&space;\in&space;[1.75,&space;1.95]&space;\text{&space;GeV},\;&space;Q^2&space;\in&space;[0.65,&space;0.7]&space;\text{&space;GeV}^2,\\&space;{&space;}\;\;\;\;\;\cos{\theta}&space;\in&space;[-1,&space;1],\;&space;\varphi&space;\in&space;[15,&space;45]&space;{&space;}^{\degree}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?E&space;=&space;6.5&space;\text{&space;GeV},&space;\;W&space;\in&space;[1.75,&space;1.95]&space;\text{&space;GeV},\;&space;Q^2&space;\in&space;[0.65,&space;0.7]&space;\text{&space;GeV}^2,\\&space;{&space;}\;\;\;\;\;\cos{\theta}&space;\in&space;[-1,&space;1],\;&space;\varphi&space;\in&space;[15,&space;45]&space;{&space;}^{\degree}" title="E = 6.5 \text{ GeV}, \;W \in [1.75, 1.95] \text{ GeV},\; Q^2 \in [0.65, 0.7] \text{ GeV}^2,\\ { }\;\;\;\;\;\cos{\theta} \in [-1, 1],\; \varphi \in [15, 45] { }^{\degree}" /></a>

## Point-like calculation
To work with this mode, when starting the program, you need to fill the boundary values with the same values for each variable:
* --W_min=X1
* --W_max=X1
* --Q2_min=X3
* --Q2_max=X3
* --cos_min=X5
* --cos_max=X5
* --phi_min=X7
* --phi_max=X7

The default values of these variables are zero, so you always need to specify the point you want.
 
> An example: ./start -e 6.5 --W_min=1.8 --W_max=1.8 --Q2_min=0.9 --Q2_max=0.9 --cos_min=0.1 --cos_max=0.1 --phi_min=180 --phi_max=180
> 
> Meaning: <a href="https://www.codecogs.com/eqnedit.php?latex=K\Lambda" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K\Lambda" title="K\Lambda" /></a> channel, point-like cross section calculation, <a href="https://www.codecogs.com/eqnedit.php?latex=E_{beam}&space;=&space;6.5\;\text{GeV},\;W&space;=&space;1.8\;\text{GeV},\;Q^2&space;=&space;0.9\;\text{GeV$^2$},\;\cos{\theta}=0.1,\;\varphi&space;=&space;180^{\degree}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?E_{beam}&space;=&space;6.5\;\text{GeV},\;W&space;=&space;1.8\;\text{GeV},\;Q^2&space;=&space;0.9\;\text{GeV$^2$},\;\cos{\theta}=0.1,\;\varphi&space;=&space;180^{\degree}" title="E_{beam} = 6.5\;\text{GeV},\;W = 1.8\;\text{GeV},\;Q^2 = 0.9\;\text{GeV$^2$},\;\cos{\theta}=0.1,\;\varphi = 180^{\degree}" /></a>
