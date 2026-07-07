"""Air<->vacuum wavelength conversion, atmospheric continuum absorption,
transmission, and airmass.

Ports the atmosphere-related routines in `src/gsetc.c` (air2vac iteration,
`gsAtmTrans`, `gsAtmContOp`, airmass from zenith angle), using the extracted
`atm_trans_kp` / `mk_trans_3mm` lookup tables. Not yet implemented.
"""
