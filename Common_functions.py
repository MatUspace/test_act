# Created:          13.07.2022
# Last Revision:    13.07.2022
# Authors:          ?, Mathieu Udriot
# Emails:           ?, mathieu.udriot@epfl.ch
# Description:      Adapted from TCAT code Common_functions (version 2022, eSpace), used to compute delta_v in orbit

from numpy.linalg import norm
from Manoeuvre import *
from poliastro.twobody import Orbit
from poliastro.core import perturbations

# global parameters for space debris index computation
ALTITUDE_ATMOSPHERE_LIMIT = 200 * u.km

def instant_orbital_velocity(orbit, radius):
    """ Returns instantaneous orbital velocity at a particular distance from attractor.
    Used during delta v calculations.

     Args:
        orbit (poliastro.twobody.Orbit): orbit
        radius (u.<distance unit>): distance to the center of the attractor

    Return:
        (u.deg): orbital speed at distance given in argument
    """
    # Check if radius is smaller than apogee plus 1m to account for rounding errors.
    if radius > orbit.r_a + 1 * u.m:
        raise Exception('Unattainable radius specified.', radius, orbit)
    # compute speed
    body = orbit.attractor
    speed = np.sqrt(body.k * (2./radius - 1/orbit.a)).to(u.m / u.s)
    return speed.to(u.m/u.s)
    
def inclination_change_delta_v(initial_orbit, final_orbit):
    """ Returns delta_v necessary to change inclination from initial to final orbital plane.
    This only includes inclination change.

    Args:
        initial_orbit (poliastro.twobody.Orbit): initial orbit
        final_orbit (poliastro.twobody.Orbit): final orbit

    Return:
        (u.m / u.s): delta v necessary for inclination change
    """
    body = initial_orbit.attractor
    delta_i = final_orbit.inc - initial_orbit.inc
    ecc = initial_orbit.ecc
    w = initial_orbit.argp
    a = max(initial_orbit.a, final_orbit.a)
    n = np.sqrt(body.k / a**3)
    f = -initial_orbit.argp
    delta_v = abs(2 * np.sin(delta_i / 2.) * (np.sqrt(1-ecc**2)*np.cos(w + f)*n*a) / (1 + ecc * np.cos(f)))
    return delta_v.to(u.m / u.s)

def high_thrust_delta_v(initial_orbit, final_orbit, initial_mass, mean_thrust, isp):
    """Returns the delta v necessary to perform an orbit change, assuming impulsive maneuvers.
    This takes into account the transfer from one elliptical orbit to another.
    This takes into account possible inclination changes, performed during the adequate impulse.
    This neglects argument of periapsis changes.
    
    Args:
        initial_orbit (poliastro.twobody.Orbit): initial orbit
        final_orbit (poliastro.twobody.Orbit): final orbit
        initial_mass (u.kg) assumed servicer mass at start of maneuver
        mean_thrust (u.N): assumed thrust at start of maneuver
        isp (u.s): assumed isp, used to estimate manoeuvre duration

    Return:
        (u.m / u.s): total delta v to reach final orbit
        (poliastro.twobody.Orbit): transfer orbit if applicable
        (u.m / u.s): first impulse
        (u.m / u.s): second impulse
        (u.day): first impulse duration
        (u.day): second impulse duration
        (u.day): total orbit change duration
    """
    manoeuvres = []
    # compute delta v for inclination change and find if inclination needs to be done during first or second impulse
    inc_delta_v = inclination_change_delta_v(initial_orbit, final_orbit)
    if initial_orbit.a > final_orbit.a:
        first_inc_delta_v = inc_delta_v
        second_inc_delta_v = 0. * u.m/u.s
    else:
        first_inc_delta_v = 0. * u.m/u.s
        second_inc_delta_v = inc_delta_v

    # let's simplify the problem by neglecting argument of periapsis changes
    # we suppose arguments of periapsis are either aligned or opposed
    # TODO: introduce argument of periapsis changes
    first_burn_radius = initial_orbit.r_a
    if abs(final_orbit.argp - initial_orbit.argp) < 180. * u.deg:
        second_burn_radius = final_orbit.r_p
    else:
        second_burn_radius = final_orbit.r_a

    # find transfer orbit, neglecting argument of periapsis change
    a = (first_burn_radius + second_burn_radius) / 2.
    ecc = abs(first_burn_radius - second_burn_radius) / (first_burn_radius + second_burn_radius)
    transfer_orbit = Orbit.from_classical(final_orbit.attractor, a, ecc, final_orbit.inc, final_orbit.raan,
                                          final_orbit.argp, final_orbit.nu, final_orbit.epoch)

    if final_orbit.attractor != initial_orbit.attractor:
        raise ValueError("Initial and final orbits have different attractors.")

    # first burn
    v_i_1 = instant_orbital_velocity(initial_orbit, first_burn_radius)
    v_f_1 = instant_orbital_velocity(transfer_orbit, first_burn_radius)
    delta_v_1 = np.sqrt((v_f_1 - v_i_1)**2 + first_inc_delta_v**2)
    manoeuvre = Manoeuvre(delta_v_1,"first high-trust dV")
    burned_mass = manoeuvre.compute_burn_duration(initial_mass, mean_thrust, isp)
    manoeuvres.append(manoeuvre)

    # second burn
    # TODO add condition on second burn if altitude is lower than ALTITUDE_ATMOSPHERIC_LIMIT
    v_i_2 = instant_orbital_velocity(transfer_orbit, second_burn_radius)
    v_f_2 = instant_orbital_velocity(final_orbit, second_burn_radius)
    delta_v_2 = np.sqrt((v_f_2 - v_i_2)**2 + second_inc_delta_v**2)
    manoeuvre = Manoeuvre(delta_v_2,"second high-trust dV")
    burned_mass = burned_mass + manoeuvre.compute_burn_duration(initial_mass, mean_thrust, isp)
    manoeuvres.append(manoeuvre)

    transfer_duration = (transfer_orbit.period / 2).to(u.day)


    return manoeuvres, transfer_duration, transfer_orbit, burned_mass


# from poliastro example https://docs.poliastro.space/en/stable/examples/Natural%20and%20artificial%20perturbations.html
# def natural_decay(t0, state, k, R, C_D, A_over_m, H0, rho0, orbit):
#     # in progress
#     """
#     Args: 
        
#         orbit: orbit in LEO on which the spacecraft is left at end of mission (either the operational orbit if no EOL startegy, or the disposal orbit if there is a disposal manoeuvre)
#     """

#     atmosphere_orbit = Orbit.from_classical(orbit.attractor, ALTITUDE_ATMOSPHERE_LIMIT+orbit.attractor.R, orbit.ecc,
#                                             orbit.inc, orbit.raan, orbit.argp, orbit.nu)

#     if orbit.ecc > 0.1:
#         raise Exception('Use of Edelbaum not valid for elliptic orbits')
        
#     # compute necessary inputs for Edelbaum formulations
#     initial_radius = (orbit.r_a + orbit.r_p) / 2
#     final_radius = (atmosphere_orbit.r_a + atmosphere_orbit.r_p) / 2
#     v_0 = instant_orbital_velocity(orbit, initial_radius)
#     v_f = instant_orbital_velocity(atmosphere_orbit, final_radius)

#     Delta_v = v_f - v_0

#     a_drag = perturbations.atmospheric_drag_exponential(t0, state, k, R, C_D, A_over_m, H0, rho0)
#     atm_perturbation = np.array([0, 0, 0, a_drag[0], a_drag[1], a_drag[2]])

#     decay_duration = (Delta_v / norm(a_drag)).to(u.year)

#     du_kep = func_twobody(t0, state, k)

#     f = du_kep + atm_perturbation
#     decay = orbit.to_ephem(strategy=EpochBounds(orbit.epoch, orbit.epoch + decay_duration.to(u.s)), method=CowellPropagator(f=f))
