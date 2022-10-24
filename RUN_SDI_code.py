# Created:          15.07.2022
# Last Revision:    15.07.2022
# Authors:          Mathieu Udriot
# Emails:           mathieu.udriot@epfl.ch
# Description:      Script for computation needed for the space debris index to be used in ACT, inputs

from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from space_debris_index import *

def main():

    print('Creating inputs...')
    # input parameters

    # TODO allow users to define more than 1 operational orbit, ideally would use more of TCAT to compute impact at each manoeuvres, also with changes in inclination ? (can be neglected for launchers but not for active spacecrafts (kickstage, ADR servicers))
    # TODO include collision avoidance manoeuvres and passivation in the assessment ?
    # TODO add check with propellant mass if maneuvres can be performed (would be solved by using TCAT directly)
    # TODO add other means of manoeuvring else than propulsive (drag sails, tumbling ?)

    mass = 1195 * u.kg
    if mass.value <= 0:
        raise ValueError("Mass muss be a positive number (kg).")

    # TODO the cross section can be computed from the CROC tool (from ESA) based on the dimensions --> what is possible to do ?
    cross_section = 15.5 * u.m ** 2
    if cross_section.value <= 0:
        raise ValueError("Cross section muss be a positive number (m^2).")
    
    starting_epoch = Time("2018-01-01 12:00:00", scale="tdb")
    op_ending_epoch = Time("2030-01-01 12:00:00", scale="tdb")
    op_duration = (op_ending_epoch - starting_epoch).to(u.year)

    mean_thrust = 8 * u.N
    Isp = 210 * u.s

    apogee_object_op = 800 * u.km
    perigee_object_op = 800 * u.km
    if perigee_object_op <= 0:
        raise ValueError("Operational perigee muss be a positive number (km).")
    elif apogee_object_op < perigee_object_op:
        raise ValueError("Operational apogee muss be a larger or equal to perigee (km).")

    inc_object_op = 98 * u.deg
    if inc_object_op >= 180 * u.deg:
        raise ValueError("Operational inclination not in the range 0 <= inc < 180.")
    elif inc_object_op < 0 * u.deg:
        raise ValueError("Operational inclination not in the range 0 <= inc < 180.")

    a_op = (apogee_object_op + perigee_object_op) / 2 + Earth.R
    ecc_op = (apogee_object_op + Earth.R - a_op) / a_op * u.one
    operational_orbit = Orbit.from_classical(Earth, 
                                            a_op, 
                                            ecc_op, 
                                            inc_object_op, 
                                            0. * u.deg,
                                            0. * u.deg,
                                            180. * u.deg,
                                            starting_epoch)
    
    EOL_manoeuvre = True

    if EOL_manoeuvre == True:
        # TODO account for Post Mission Disposal (PMD) Success Rate by computing the impacts with and without successful manoeuvres and sum the two
        PMD_success = 0.9
        
        apogee_object_disp = 800 * u.km
        perigee_object_disp = 600 * u.km
        if perigee_object_disp <= 0:
            raise ValueError("Disposal perigee muss be a positive number (km).")
        elif apogee_object_disp < perigee_object_disp:
            raise ValueError("Disposal apogee muss be a larger or equal to perigee (km).")

        inc_object_disp = 98 * u.deg
        if inc_object_disp >= 180 * u.deg:
            raise ValueError("Disposal inclination not in the range 0 <= inc < 180.")
        elif inc_object_disp < 0 * u.deg:
            raise ValueError("Disposal inclination not in the range 0 <= inc < 180.")

        a_disp = (apogee_object_disp + perigee_object_disp) / 2 + Earth.R
        ecc_disp = (apogee_object_disp + Earth.R - a_disp) / a_disp * u.one
        disposal_orbit = Orbit.from_classical(Earth, 
                                                a_disp, 
                                                ecc_disp, 
                                                inc_object_disp, 
                                                0. * u.deg,
                                                0. * u.deg,
                                                180. * u.deg,
                                                starting_epoch)
    
    # case without manoeuvre
    else:
        PMD_success = 0
        disposal_orbit = operational_orbit
    
    # input .csv file for characterization factor
    CF_file = np.genfromtxt(f'space_debris_CF_for_code.csv', delimiter=",", skip_header=1)
    # input .csv file for natural decay
    reduced_lifetime_file = np.genfromtxt(f'reduced_lifetime.csv', delimiter=",", skip_header=2)

    SDI_results = SDI_compute(CF_file, reduced_lifetime_file, mass, cross_section, op_duration, mean_thrust, Isp, operational_orbit, EOL_manoeuvre, disposal_orbit)

def SDI_compute(CF_file, reduced_lifetime_file, mass, cross_section, op_duration, mean_thrust, Isp, operational_orbit, EOL_manoeuvre, disposal_orbit):
    # Impact score computation
    print('Start finding the orbital case and computing impact score...')

    # #1 Case if operational perigee is lower than the atmosphere limit, no natural decay hereafter
    if (operational_orbit.r_p - Earth.R) < ALTITUDE_ATMOSPHERE_LIMIT:
        # Case if object is already completely in the atmosphere
        if (operational_orbit.r_a - Earth.R) < ALTITUDE_ATMOSPHERE_LIMIT:
            print("Object will reenter directly.")
            impact_score = 0 * u.year *u.pot_fragments
            op_impact_percentage = 0
            disp_maneuver_impact_percentage = 0
            natural_impact_percentage = 0
            print("Computed space debris impact score is", impact_score, ".")

            results = {"Space_Debris_Index": impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
            return results
        else:
            # Case if apogee is in space and perigee in the atmosphere: direct reentry
            print("Direct reentry.")
            # Compute impact of reentry: decompose elliptical orbit and use time_to_anomaly after finding LTAN from position meaning altitude
            disposal_impact = elliptical_orbit_decomposition(CF_file, operational_orbit, mass)

            impact_score = disposal_impact*mass*cross_section
            op_impact_percentage = 0
            disp_maneuver_impact_percentage = 100
            natural_impact_percentage = 0
            print("Computed space debris impact score is", "{:.3f}".format(impact_score), ". 0 percent operational impact, 100 percent disposal impact.")

            results = {"Space_Debris_Index": impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
            return results

    # #2 Case if object's perigee is higher than LEO
    elif (operational_orbit.r_p - Earth.R) > ALTITUDE_LEO_LIMIT:
        print("Operational orbit is higher than LEO, no debris impact computed for the operational phase.")
        CF_op = 0
        op_impact_percentage = 0
        # Check disposal orbit
        if (disposal_orbit.r_p - Earth.R) > ALTITUDE_LEO_LIMIT:
            print("Graveyard orbit outside of LEO, no debris impact computed for the disposal phase.")
            # TODO add decay from above LEO limit ?
            impact_score = 0 * u.year *u.pot_fragments
            disp_maneuver_impact_percentage = 0
            natural_impact_percentage = 0
            print("Computed space debris impact score is", impact_score, ".")

            results = {"Space_Debris_Index": impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
            return results
        else:
            print("Reentry manoeuvre from operational orbit higher than LEO.")
            manoeuvres, transfer_duration, transfer_orbit, burned_mass = high_thrust_delta_v(operational_orbit, disposal_orbit, mass, mean_thrust, Isp)
            # Compute impact of disposal manoeuvre: decompose elliptical orbit and use time_to_anomaly after finding LTAN from position meaning altitude
            print("Disposal manoeuvre")
            disposal_impact = elliptical_orbit_decomposition(CF_file, transfer_orbit, mass - burned_mass)*(mass - burned_mass)*cross_section

            # TODO switch with natural decay w/o disposal, weighted by PMD success rate of disposal maneuvre ?
            print("Natural decay")
            natural_decay_time, natural_decay_impact = natural_decay(reduced_lifetime_file, CF_file, disposal_orbit, cross_section, mass, transfer_duration, op_duration)
            
            total_impact_score = disposal_impact + natural_decay_impact

            natural_impact_percentage = natural_decay_impact/total_impact_score*100
            disp_maneuver_impact_percentage = 100 - natural_impact_percentage
            print("Computed space debris impact score is", "{:.3f}".format(total_impact_score), ". 0 percent operational impact, 100 percent disposal impact. Of which", "{:.3f}".format(disp_maneuver_impact_percentage), "percent from the disposal manoeuvre", "{:.3f}".format(natural_impact_percentage), "percent from the natural decay.")
            
            results = {"Space_Debris_Index": total_impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
            return results

    # #3 Cases if object has perigee in LEO
    # Case of LEO circular orbit
    elif operational_orbit.r_a == operational_orbit.r_p:
        print("LEO circular")
        CF_op = get_characterization_factor(CF_file, (operational_orbit.r_p - Earth.R), operational_orbit.inc)
        OP_impact = cross_section*mass*CF_op*alpha_param(mass)*op_duration

        # TODO switch with natural decay w/o disposal, weighted by PMD success rate of disposal maneuvre ?
        # disposal manoeuvre
        if EOL_manoeuvre:
            print("Disposal manoeuvre")
            manoeuvres, transfer_duration, transfer_orbit, burned_mass = high_thrust_delta_v(operational_orbit, disposal_orbit, mass, mean_thrust, Isp)
            # Compute impact of disposal manoeuvre
            disposal_impact = elliptical_orbit_decomposition(CF_file, transfer_orbit, mass - burned_mass) # intermediate impact
        else:
            disposal_impact = 0 * u.pot_fragments * u.year * u.kg **(-1) * u.m **(-2) # intermediate impact
            transfer_duration = 0 * u.day
            burned_mass = 0 * u.kg

        # natural decay impact
        print("Natural decay")
        natural_decay_time, natural_decay_impact = natural_decay(reduced_lifetime_file, CF_file, disposal_orbit, cross_section, mass, transfer_duration, op_duration)

        total_impact_score = OP_impact + cross_section*(mass - burned_mass)*disposal_impact + natural_decay_impact

        natural_impact_percentage = natural_decay_impact/total_impact_score*100
        op_impact_percentage = OP_impact/total_impact_score*100
        disp_maneuver_impact_percentage = 100 - natural_impact_percentage - op_impact_percentage

        print("Computed space debris impact score is", "{:.3f}".format(total_impact_score), ".", "{:.3f}".format(op_impact_percentage), "percent from operations", "{:.3f}".format(disp_maneuver_impact_percentage), "percent from disposal manoeuvre", "{:.3f}".format(natural_impact_percentage), "percent from natural decay.")
        
        results = {"Space_Debris_Index": total_impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
        return results

    # Cases of elliptical operational orbit
    else:
        print("Elliptical operational orbit partially higher than LEO or completely in LEO.")

        # Compute impact of half elliptical operational orbit, times two for the complete impact of one pass in LEO, times the number of orbits during the operational lifetime
        operational_impact = cross_section*mass*2*elliptical_orbit_decomposition(CF_file, operational_orbit, mass)*op_duration/operational_orbit.period.to(u.year)

        # TODO switch with natural decay w/o disposal, weighted by PMD success rate of disposal maneuvre ?
        # disposal manoeuvre
        if EOL_manoeuvre:
            print("Disposal manoeuvre")
            manoeuvres, transfer_duration, transfer_orbit, burned_mass = high_thrust_delta_v(operational_orbit, disposal_orbit, mass, mean_thrust, Isp)
            # decompose elliptical disposal  orbit and use time_to_anomaly after finding LTAN from position meaning altitude
            disposal_impact = cross_section*(mass - burned_mass)*elliptical_orbit_decomposition(CF_file, transfer_orbit, mass - burned_mass) # intermediate impact
        else:
            disposal_impact = 0 * u.pot_fragments * u.year * u.kg **(-1) * u.m **(-2) # intermediate impact
            transfer_duration = 0 * u.day
            burned_mass = 0 * u.kg
        
        print("Natural decay")
        natural_decay_time, natural_decay_impact = natural_decay(reduced_lifetime_file, CF_file, disposal_orbit, cross_section, mass, transfer_duration, op_duration)

        total_impact_score = operational_impact + disposal_impact + natural_decay_impact

        natural_impact_percentage = natural_decay_impact/total_impact_score*100
        op_impact_percentage = operational_impact/total_impact_score*100
        disp_maneuver_impact_percentage = 100 - natural_impact_percentage - op_impact_percentage

        print("Computed space debris impact score is", "{:.3f}".format(total_impact_score), ".", "{:.3f}".format(op_impact_percentage), "percent from operations", "{:.3f}".format(disp_maneuver_impact_percentage), "percent from disposal manoeuvre", "{:.3f}".format(natural_impact_percentage), "percent from natural decay.")
        
        results = {"Space_Debris_Index": total_impact_score, "Operational_percentage": op_impact_percentage, "Disposal_manoeuvre_percentage": disp_maneuver_impact_percentage, "Natural_decay_percentage": natural_impact_percentage}
        return results

if __name__ == "__main__":
    main()