# Created:          27.09.2022
# Last Revision:    27.09.2022
# Authors:          Mathieu Udriot
# Emails:           mathieu.udriot@epfl.ch
# Description:      Script to read input .csv file describing the trajectory and thrust curve needed to compute the atmoshperic impacts of a launch vehicle, to be used in ACT

from scipy import interpolate
import numpy as np
import astropy.units as u
from matplotlib import pyplot

from Common_functions import *
# atmospheric limit is defined in common functions, limit is set to split impacts on the atmosphere and impacts on the space environment
# atmospheric layers are defined below:
ATM_LIM_TROPOSPHERE = 10 * u.km
ATM_LIM_OZONE_LOW = 20 * u.km
ATM_LIM_OZONE_HIGH = 30 * u.km
ATM_LIM_STRATOSPHERE = 50 * u.km
ATM_LIM_MESOSPHERE = 85 * u.km

PATH_CSV__THRUST_CURVES = "atm_thrust_curves/"
PATH_CSV_TRAJECTORIES = "atm_trajectories/"

def read_inputs():
    
    launcher = "Ariane_5"
    engine = "Vulcain"

    # Trajectory up to atmospheric limit, will not be used above
    raw_trajectory = np.genfromtxt(f'{PATH_CSV_TRAJECTORIES}input_traj_{launcher}.csv', delimiter=",", skip_header=2)

    # thrust curve, code should be run with the main engine, or the booster, or the secondary engine, when needed
    raw_thrust_curve = np.genfromtxt(f'{PATH_CSV__THRUST_CURVES}thrust_curve_{engine}.csv', delimiter=",", skip_header=2)

    # ignition and cutoff altitude
    

    # if several engines used in parallel (Prometheus), this scales the emissions too. From input.
    engine_s = 1

    prop_type = 'LH2/LOx'

    return
