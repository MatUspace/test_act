# Created:          13.07.2022
# Last Revision:    13.07.2022
# Authors:          ?, Mathieu Udriot
# Emails:           ?, mathieu.udriot@epfl.ch
# Description:      Adapted from TCAT code Manoeuvre.py (version 2022, eSpace)

import numpy as np
from astropy import units as u
from astropy import constants as const

class Manoeuvre:
    """ Class representing a manoeuvre. This is used to simplify computations of thrust, mass and durations.

    Args:
        delta_v (u.<speed unit>): delta v for the manoeuvre

    Attributes:
        delta_v (u.<speed unit>): delta v for the manoeuvre
        burn_duration (u.<time unit>): duration of the burn for the manoeuvre
    """
    def __init__(self, delta_v, id = ""):
        self.delta_v = delta_v
        self.burn_duration = 0. * u.minute
        self.id = id

    def get_delta_v(self):
        return self.delta_v

    def compute_burn_duration(self, initial_mass, mean_thrust, isp):
        final_mass = initial_mass / np.exp((self.delta_v.to(u.meter / u.second) / const.g0 / isp.to(u.second)).value)
        mean_mass = (final_mass + initial_mass) / 2
        self.burn_duration = mean_mass / mean_thrust * self.delta_v
        self.burn_duration = self.burn_duration.to(u.s)
        return initial_mass - final_mass

    def get_burn_duration(self, duty_cycle=1.):
        return self.burn_duration / duty_cycle

    def convert_time_for_print(dt):
        """ Returns the time in the lowest unit with decimal > 0.
        Typ:    dt = 135s --> returns: 2.25 minute
                dt = 3600s --> returns: 1 hour
        Args:
            dt (u.sec, u.minute, u.day, etc...): duration to covert in a printable unit
        """
        if dt > 1. * u.year:
            dt = dt.to(u.year)
        elif dt > 1. * u.day:
            dt = dt.to(u.day)
        elif dt > 1 * u.hour:
            dt = dt.to(u.hour)
        elif dt > 1 * u.minute:
            dt = dt.to(u.minute)
        elif dt > 1 * u.s:
            dt = dt.to(u.s)
        else:
            dt = dt.to(u.ms)
        return dt

    def __str__(self):
        duration_print = self.convert_time_for_print(self.burn_duration)
        return (f"\u0394V: {self.delta_v.to(u.m/u.s):.1f}, \u0394t {duration_print:.1f}, {self.id}")