import sys
sys.path.append('../src')

from ms_toolbox import tools as tls

import matplotlib.pyplot as plt
import numpy as np

from ms_toolbox import gfd
from ms_toolbox import rgmodel as rgm
from pycurrents.num import spectra
from pycurrents.file.matfile import loadmatbunch

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp2d as int2d

from netCDF4 import Dataset
from pycurrents.system import Bunch
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import RegularGridInterpolator as RGI

import cmocean
import matplotlib.cm as cm
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'

###### Function to load in variables ######

figpath = '../small_islands/compare/figures/'

fcor = gfd.coriolis(8)

def load_var(varpath,file_name1,gp):
    a = np.load(varpath + file_name1)
    b = Bunch()
    b.cp_m = a["cp"]
    b.cpf = a["cpf"]
    b.cp_t = np.sqrt(gp * 100)  #LH00
    b.epsilon_L = a["eps"]

    return b

fpath1 = '../small_islands/compare/data_subsection/'
file_name1 = "azimuthal_cp_rad_f.npz"
################# Define Paths ###################

anpath = fpath1 + 'r50/'
####################################################

rg002 = load_var(anpath,file_name1,.09)

################# Define Paths ###################

anpath = fpath1 + 'r75/'
####################################################

rg009 = load_var(anpath,file_name1,.09)


################# Define Paths ###################

anpath = fpath1 + 'r100/'
####################################################

rg013 = load_var(anpath,file_name1,.09)


################# Define Paths ###################

anpath = fpath1 + 'r150/'
####################################################

rg017 = load_var(anpath,file_name1,.09)

################# Define Paths ###################

anpath = fpath1 + 'r200/'
####################################################

rg025 = load_var(anpath,file_name1,.09)

####################################################


file_name2 = "analytical/theory_tw.npz"
theo =  np.load(fpath1 + file_name2)
eps_t = theo["eps"]
sigf_t = theo["sigf"]


cp_m = (rg002.cp_m,rg009.cp_m,rg013.cp_m,rg017.cp_m,rg025.cp_m)
epsilon_L = (rg002.epsilon_L,rg009.epsilon_L,rg013.epsilon_L,rg017.epsilon_L,rg025.epsilon_L)

cpf = (rg002.cpf,rg009.cpf,rg013.cpf,rg017.cpf,rg025.cpf)
fcor = gfd.coriolis(8)
a = 50e3
af = fcor[0]*a

fig,ax = plt.subplots(figsize=(10,5))
ax.plot(epsilon_L,np.abs(cpf),'.',color="deepskyblue",zorder=10,markersize=15)
ax.plot(eps_t,sigf_t,color="blue",zorder=5)

ax.set_xlim((.07,2))
ax.set_ylim((.5,1))
ax.set_xlabel("$\\epsilon_L$")
ax.set_xlabel("$\\epsilon_L$",fontsize=20)
ax.set_ylabel("Angular frequency / $f$ [$\\frac{\\partial \\theta}{\\partial t} / f$]",fontsize=20)
ax.tick_params(axis="both",labelsize=18)




fig.savefig(figpath + 'FIGURE_11' + ".png", dpi=300, bbox_inches='tight',orientation="horizontal")
