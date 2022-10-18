# Test_ACT
## Space debris Index (SDI) code

Run the space debris index (SDI) code by using the command: _python RUN_SDI_code.py SDI_inputs.json_

For now this inputs have to be changed in RUN_SDI_code.py but soon from the SDI_inputs.json file

The code starts with collecting the inputs and then uses three cases with switches:
  1. The operational perigee is below the atmospheric limit (set to 200km) which means the object is considered re-entered.
    o If the apogee is also lower than 200km, the object is completely in the atmosphere so there is zero impact.
    o Else, it’s a direct reentry (only impacts from the reentry).
  2. The operational orbit perigee is above the limit of the protected LEO region (2000km), meaning no impact is computed during operations because there are no CFs higher than 2000km.
    o If the perigee of the disposal orbit (target of the EOLM) is above 2000km (case of a graveyard orbit), no impact is computed for the disposal phase either.
    o In the case of a reentry manoeuvre, the impact of the manoeuvre and the natural decay from the disposal orbit, down to the atmospheric limit is computed.
  3. The object has its perigee in LEO:
    o Case of a circular operational orbit
    o Case of an elliptical operational orbit

The _main_ method then calls methods from space_debris_index.py and common methods for orbital mechanics from common_functions.py and manoeuvre.py (both defined in TCAT)
The code uses two .csv files for its calculation. 

## Atmospheric emissions of launchers
All files and fodlers starting with atm_ are in progress for the calculation of the atmospheric emissions during a launch.
