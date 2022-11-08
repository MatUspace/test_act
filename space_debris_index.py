# Created:          12.07.2022
# Last Revision:    12.07.2022
# Authors:          Mathieu Udriot
# Emails:           mathieu.udriot@epfl.ch
# Description:      Script for computation needed for the space debris index to be used in ACT

from astropy import units as u
from poliastro.bodies import Earth
from Common_functions import *


# New astropy unit
u.pot_fragments = u.def_unit("potential_fragments")

# global parameters for space debris index computation
ALTITUDE_LEO_LIMIT = 2000 * u.km
ALTITUDE_INCREMENT = 50 * u.km
INCLINATION_INCREMENT = 2 * u.deg
TIME_INTERVAL_LIMIT = 200 * u.year
RESIDUAL_TIME_IADC_GUIDELINE = 25 * u.year
RESIDUAL_TIME_FCC_GUIDELINE = 5 * u.year

PRINT_BOOL = False

""" A space object can be a spacecraft (servicers, upper stages, kick stages, satellites (active or defunct)) or any type of debris

The space debris index has been created as a LCA indicator (T. Maury et alli, "Space debris through the prism  of environemental performence of space systems: the case of Sentinel-3 redesigned mission", (2019))
For now the characterization factors are only available in LEO between 200 and 2000 km
Orbital cells are defined by altitude and inclination

"""

# methods
def alpha_param(mass):
    """ To find the alpha coefficient depending on the mass of the object 
    T. Maury et alli, "Space debris through the prism of the environemntal performance of space systems: the case of Sentinel-3 redesigned mission", 2019
    
    Values from the NASA break-up model
    """

    if mass.value <= 500:
        alpha = 1.3 * u.kg **(-1)
    elif mass.value <= 1000:
        alpha = 1 * u.kg **(-1)
    elif mass.value <= 1500:
        alpha = 0.86 * u.kg **(-1)
    elif mass.value <= 2000:
        alpha = 0.8 * u.kg **(-1)
    else:
        alpha = 0.77 * u.kg **(-1)

    return alpha

def natural_decay(reduced_lifetime_file, CF_file, initial_orbit, cross_section, mass, disposal_time, op_time, type):
    """
    To compute natural decay from a LEO orbit down to the ALTITUDE_ATMOSPHERE_LIMIT

    Orbital stages left above LEO (perigee > ALTITUDE_LEO_LIMIT) are considered to stay there for too long to be assessed (in fact there are no CFs for now so their impact score would anyway be 0)

    cells in the reduced lifetime have Delta altitude equal to ALTITUDE_INCREMENT 

    :param initial_orbit: Orbit on which the spacecraft starts decaying
    :type initial_orbit: poliastro.twobody.Orbit
    :param cross_section: Randomly tumbling cross section
    :type cross_section: u*m**2
    :param mass: Mass of the spacecraft
    :type mass: u*kg
    :param disposal_time: Transfer time of the disposal manoeuvre
    :type disposal_time: u*day
    :param op_time: Duration of the operational phase
    :type op_time: u*year
    :param type: Type of impact: after successful or unsuccessful EOLM
    :type type: boolean (True or False)

    Output:
    :cumulated_time [yrs]
    :cumulated natural decay impact [pot. fragments * yrs]
    """
    perigee = initial_orbit.r_p - Earth.R
    ecc = initial_orbit.ecc
    inc = initial_orbit.inc

    A_over_m = cross_section / mass
    decay_orbit_impact = 0

    # time interval is caped at 200 yrs because the characterization factors are computed using this limit. So the cumulative time (including operations and disposal manoeuvre) shall not go over 200 yrs (TIME_INTERVAL_LIMIT).
    cumulated_time = op_time + disposal_time.to(u.year)
    total_disposal_time = disposal_time.to(u.year)
    
    # TODO add decay from above 2000km in case it enters the LEO protected region and add if statement to start impact integration only then ? Might not make sense due to really slow
    # find perigee index (in this case between 0 (200 km) and 36 (2000km), with increment of 50km)
    index_peri = int((round(perigee.value*2)/2 - ALTITUDE_ATMOSPHERE_LIMIT.value)/ALTITUDE_INCREMENT.value)

    if index_peri < 0:
        print("Direct reentry, no impact from natural decay.")
        decay_orbit_impact = 0 * u.year * u.pot_fragments
        return cumulated_time, decay_orbit_impact

    # find eccentricity index
    if ecc >= 0.8:
        index_ecc = 7
    elif ecc >= 0.4:
        index_ecc = 6
    elif ecc >= 0.2:
        index_ecc = 5
    elif ecc >= 0.1:
        index_ecc = 4
    elif ecc >= 0.05:
        index_ecc = 3
    elif ecc >= 0.02:
        index_ecc = 2
    elif ecc >= 0.01:
        index_ecc = 1
    else:
        index_ecc = 0

    # case of a circular orbit, assumed to stay circular, decay of time
    if ecc <0.01:
        print("natural decay from circular orbit")
        while index_peri > 0: # so it doesn't read first column
            
            cell_time = (reduced_lifetime_file[0, index_peri] - reduced_lifetime_file[0, index_peri - 1])/A_over_m.value * u.year

            total_disposal_time += cell_time
            cumulated_time += cell_time

            if cumulated_time > TIME_INTERVAL_LIMIT:
                decay_orbit_impact += get_characterization_factor(CF_file, perigee, inc)*(TIME_INTERVAL_LIMIT - (cumulated_time - cell_time))*alpha_param(mass)*cross_section*mass
                print("/!\ Time interval limit reached !")
                break
            else:
                decay_orbit_impact += get_characterization_factor(CF_file, perigee, inc)*cell_time*alpha_param(mass)*cross_section*mass

            perigee = perigee - ALTITUDE_INCREMENT
            index_peri -= 1
    else:
        # first big assumption that only the eccentricity is decreasing because most braking is done at perigee which lowers apogee, until the orbit is circular. 
        # But it is decreased only after a given time cell_time in the same orbit after which it jumps to the next eccentricity value.
        # computing a new eccentricity (so apogee) after each passage at perigee would be computing intensive I think
        print("natural decay from elliptical orbit")
        while index_ecc > 0:
            cell_time = (reduced_lifetime_file[index_ecc, index_peri] - reduced_lifetime_file[index_ecc - 1, index_peri])/A_over_m.value * u.year
            
            total_disposal_time += cell_time
            cumulated_time += cell_time

            a = (perigee + Earth.R)*(1+ecc)/(1-ecc) #semi-major axis of ellipse

            elliptical_orbit = Orbit.from_classical(Earth, 
                                            a, 
                                            ecc, 
                                            inc, 
                                            0. * u.deg,
                                            0. * u.deg,
                                            180. * u.deg)

            if cumulated_time > TIME_INTERVAL_LIMIT:
                decay_orbit_impact += elliptical_orbit_decomposition(CF_file, elliptical_orbit, mass)*cross_section*mass*2*(TIME_INTERVAL_LIMIT - (cumulated_time - cell_time))/elliptical_orbit.period.to(u.year)
                print("/!\ Time interval limit reached !")
                break
            else:
                decay_orbit_impact += elliptical_orbit_decomposition(CF_file, elliptical_orbit, mass)*cross_section*mass*2*cell_time/elliptical_orbit.period.to(u.year)
            
            index_ecc -= 1
            ecc = reduced_lifetime_file[index_ecc, 0] * u.one

        # once orbit is circular, the decay is assumed to take place with always quasi circular orbits
        print("natural decay continues as circular orbit")
        while index_peri > 0: # so it doesn't read first column
            cell_time = (reduced_lifetime_file[0, index_peri] - reduced_lifetime_file[0, index_peri - 1])/A_over_m.value * u.year
            total_disposal_time += cell_time
            cumulated_time += cell_time

            if cumulated_time > TIME_INTERVAL_LIMIT:
                decay_orbit_impact += get_characterization_factor(CF_file, perigee, inc)*(TIME_INTERVAL_LIMIT - cumulated_time - cell_time)*alpha_param(mass)*cross_section*mass
                print("/!\ Time interval limit reached !")
                break
            else:
                decay_orbit_impact += get_characterization_factor(CF_file, perigee, inc)*cell_time*alpha_param(mass)*cross_section*mass

            perigee = perigee - ALTITUDE_INCREMENT
            index_peri -= 1
    
    if type:
        print("Disposal and residual orbital lifetime after successful EOLM equals", "{:.3f}".format(total_disposal_time))
        # flag non compliance with guidelines
        if total_disposal_time > RESIDUAL_TIME_IADC_GUIDELINE:
            print("/!\ Disposal and residual orbital lifetime after successful EOLM exceeds the IADC guideline of", RESIDUAL_TIME_IADC_GUIDELINE, "years.")
        if total_disposal_time > RESIDUAL_TIME_FCC_GUIDELINE:
            print("/!\ Disposal and residual orbital lifetime after successful EOLM exceeds the FCC guideline of", RESIDUAL_TIME_FCC_GUIDELINE, "years.")
    else:
        print("Disposal and residual orbital lifetime due to unsuccessful EOLM equals", "{:.3f}".format(total_disposal_time))
        # flag non compliance with guidelines
        if total_disposal_time > RESIDUAL_TIME_IADC_GUIDELINE:
            print("/!\ Disposal and residual orbital lifetime due to unsuccessful EOLM exceeds the IADC guideline of", RESIDUAL_TIME_IADC_GUIDELINE, "years.")
        if total_disposal_time > RESIDUAL_TIME_FCC_GUIDELINE:
            print("/!\ Disposal and residual orbital lifetime due to unsuccessful EOLM exceeds the FCC guideline of", RESIDUAL_TIME_FCC_GUIDELINE, "years.")

    return cumulated_time, decay_orbit_impact

def elliptical_orbit_decomposition(CF_file, transfer_orbit, mass):
    """ To decompose the trajectory in time spent in different orbital cell (altitude and inclination)

    Depends on global parameters defined by the characterization factors available

    /!\ For now only computes score when transfer orbit comes from higher to lower altitude
    TODO add case when transfer orbit is to bring object higher (eg. graveyard orbit > 2000km from within LEO)
    If switch at the beginning to assign global constant to condition values ? Delta altitude = +/- ALTITUDE_INCREMENT, integration boundary = ALTITUDE_ATMOSPHERE_LIMIT or ALTITUDE_LEO_LIMIT...

    return
        (yr*potential_fragments/m**2 /kg) intermediate impact score
    """

    apogee_transfer = (transfer_orbit.r_a - Earth.R)
    perigee_transfer = (transfer_orbit.r_p - Earth.R)

    if apogee_transfer > ALTITUDE_LEO_LIMIT + ALTITUDE_INCREMENT:
        altitude_temp = ALTITUDE_LEO_LIMIT 
        nu_1 = get_nu(transfer_orbit.a, transfer_orbit.ecc, altitude_temp + Earth.R)
    else:
        altitude_temp = apogee_transfer
        nu_1 = 179*np.pi/180 * u.rad

    nu_temp = nu_1

    t_temp = Orbit.time_to_anomaly(transfer_orbit, nu_temp)
    impact_score_temp = 0 * u.year* u.m **(-2)*u.pot_fragments

    # integrates impact from apogee to atmospheric limit or perigee (so during maximum half an ellipse)
    while altitude_temp - ALTITUDE_INCREMENT >= max(ALTITUDE_ATMOSPHERE_LIMIT, perigee_transfer):
        nu_next = get_nu(transfer_orbit.a, transfer_orbit.ecc, altitude_temp - ALTITUDE_INCREMENT + Earth.R)
        t_next = Orbit.time_to_anomaly(transfer_orbit, nu_next)
        Delta_t = t_temp - t_next
        Delta_t = Delta_t.to(u.year)
        impact_score_temp = impact_score_temp + Delta_t*get_characterization_factor(CF_file, altitude_temp, transfer_orbit.inc)

        if PRINT_BOOL:
            print("altitude", altitude_temp, ", anomaly nu", nu_temp, ", altitude_next", altitude_temp - ALTITUDE_INCREMENT, ", nu next", nu_next, ", Delta t", Delta_t.to(u.s))

        nu_temp = nu_next
        t_temp = t_next
        altitude_temp = altitude_temp - ALTITUDE_INCREMENT

    # Last cell, Delta_t time might be smaller
    nu_next = get_nu(transfer_orbit.a, transfer_orbit.ecc, max(ALTITUDE_ATMOSPHERE_LIMIT, perigee_transfer) + Earth.R)
    t_next = Orbit.time_to_anomaly(transfer_orbit, nu_next)
    Delta_t = t_temp - t_next
    Delta_t = Delta_t.to(u.year)
    if PRINT_BOOL:
        print("final altitude", altitude_temp, ", anomaly nu", nu_temp, ", altitude_next", max(ALTITUDE_ATMOSPHERE_LIMIT, perigee_transfer), ", nu next", nu_next, ", Delta t", Delta_t.to(u.s))

    impact_score_temp = impact_score_temp + Delta_t*get_characterization_factor(CF_file, altitude_temp, transfer_orbit.inc)

    # mass assumed not to change so same alpha parameter applies to all contributions
    return impact_score_temp*alpha_param(mass)

def get_nu(a, ecc, r):
    """ Compute the true anomaly value for a given position on the orbit
    Assumptions that only high thrust burn are performed so the transfer orbits are only half orbits
    Meaning 0 <= nu <= pi --> -1 <= cos(nu) <= 1

    Args
    a: semi-major axis of the transfer orbit (u.km)
    ecc: eccentricity of the transfer orbit (u.one)
    r: distance between the focus of attraction and the orbiting object (u.km)

    return
    (u.rad) nu: true anomaly of the orbiting object (angle defining the position of the orbiting body along the ellipse)
    """
    nu = np.arccos(round((1/ecc.value)*((a.value*(1 - ecc.value**2)/ r.value) - 1), 8)) * u.rad

    return nu

def get_characterization_factor(CF_file, altitude, inclination):
    """ Get CF value from a database
    Database has cells defined by Delta in altitude and in inclination, global parameters must be similar to accommodate those cells
    Need to round number to find the correct value in the csv file corresponding to the orbital cell

    Args:
        CF_file: .csv file imported with the characterization factors for each orbital cell
        altitude: the altitude (of the current position) on the orbit to find the index
        inclination: the inclination of the orbit to find the index

    Return:
        (Potential fragments per m^2) The charcterization factor (CF) corresponding to the input orbital cell. Potential fragments are computed by number of debris*number of fragmentation events (ref. T.Maury et alli)
    """

    # find inclination index (in this case between 1 (0 deg) and 90 (178 deg), with increment of 2deg)
    index_inc = int(np.floor(inclination.value/INCLINATION_INCREMENT.value))
    # find altitude index (in this case between 1 (200 km) and 37 (2000km), with increment of 50km)
    index_alt = int((round(altitude.value*2)/2 - ALTITUDE_ATMOSPHERE_LIMIT.value)/ALTITUDE_INCREMENT.value)

    # .csv ignores first row (header) but there is a first column with the inclination values
    characterization_factor = CF_file[index_inc, index_alt + 1] * u.m **(-2)*u.pot_fragments
    
    if PRINT_BOOL:
        print("CF:", characterization_factor)

    return characterization_factor
