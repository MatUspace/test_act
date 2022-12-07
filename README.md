# Test_ACT
## Space debris Index (SDI) code
All files and folders starting with sdi_ are for the calculation of the space debris index

Run the space debris index (SDI) code by using the command: _python sdi_run_code.py_

For now the inputs have to be changed in sdi_run_code.py. The script will be accessed with an API, from ACT, entering method sdi_main() which verifies the inputs and called SDI_compute() for the different systems.

SDI_compute() uses three cases with switches:
  1. The operational perigee is below the atmospheric limit (set to 200km) which means the object is considered re-entered.
    - If the apogee is also lower than 200km, the object is completely in the atmosphere so there is zero impact.
    - Else, it’s a direct reentry (only impacts from the reentry).
  2. The operational orbit perigee is above the limit of the protected LEO region (2000km), meaning no impact is computed during operations because there are no CFs higher than 2000km.
    - If the perigee of the disposal orbit (target of the EOLM) is above 2000km (case of a graveyard orbit), no impact is computed for the disposal phase either.
    - In the case of a reentry manoeuvre, the impact of the manoeuvre and the natural decay from the disposal orbit, down to the atmospheric limit is computed.
  3. The object has its perigee in LEO:
    - Case of a circular operational orbit
    - Case of an elliptical operational orbit

The _sdi_main()_ method calls the SDI_compute() method:
- for the orbital stage of the launch vehicle
- and if the scenario includes an ADR stage: 
    - for the ADR servicer from insertion to the target debris orbit
    - for the debris if it is not removed
    - for the servicer + debris system for the deorbitation manoeuvre

In the SDI_compute, methods from space_debris_index.py and common methods for orbital mechanics from common_functions.py and manoeuvre.py (both defined in TCAT)
The code uses two .csv files for its calculation:
- SDI characterization factors for impact calculation
- Reduced lifetime in orbit for natural decay approximations 

For each SDI computation, the results are stored and the final index is computed:
- for an LV without ADR stage: with just the SDI of the orbital stage
- for an LV with an ADR stage: by summing the SDI of the orbital stage, the ADR stage, and the servicer+debris deorbitation, and subtracting the resiudal SDI of the debris if it was not removed

The output is scaled by the number of launches.

## Atmospheric emissions of launchers
All files and folders starting with atm_ are for the calculation of the atmospheric emissions during a launch.

Run the atmospheric emissions (ATM) code by using the command: _python atm_run_code.py_

For now the inputs have to be changed in atm_run_code.py. The script will be accessed with an API, from ACT, entering method atm_main() which verifies the inputs and perform the calculations.

The _atm_main()_ method starts by interpolating the provided trajectory and the mass flow curve (extracted from the input thrust curve). They can be plotted if the Boolean "PLOTTING" is True.

The table atm_emissions_per_propellant.csv is imported, it holds values of emissions (in kg for 12 species: "CO", "CO2", "H2O", "H", "O", "OH", "N2", "NO", "Al", "HCl", "Cl", "soot (BC)") per kig of burnt propellant. Propellants in the table are:
1. LOx/RP1
2. LOx/LH2
3. LOx/LCH4
4. NTO/UDMH
5. APCP
A list of layer classes to define the atmosphere made of several layers, is created, based on gloabl parameters defining the limits.

From the interpolated trajectory, timestamps at which the engine crosses the limit between two layers are found. The code accepts trajectories going down again (for propulsive reusable LVs or for ascent that have a coast phase like Ariane 5). Between these timestamps, the mass flow curvbe is integrated to find the mass of propellant burnt in each layer. This value is used with the atm_emissions_per_propellant table to find the mass of emissions.

The output is scaled by the number of engines and launches.