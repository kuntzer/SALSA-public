# Resolution of the sky grid
resx = 40
resy = 20

# Detector characteristics
dlambda_over_lambda = 0.8
radius_psf = 15.5 # pixels
radius_tele = 0.15 # m

# Unit: V magnitude (used for clipping value, important as if mag_max is too high, 
#	the pipeline does not converge)
magnitude_min = 6
magnitude_max = 13

# Unit: ph/(px s)
sl_min = 1e-8
sl_max = 1e0

# Mirrors efficiency (up to 1 --> 70% = 0.7)
mirror_efficiency = 0.7

# Post-treatment stray light reduction (up to 1 --> 99.6% = 0.996)
SL_post_treat_reduction = 0.996

# ppm above the image is considered as contaminated by stray light
ppm_threshold = 70.

# Name of map file without markings. Ie for orbit_000.dat -> orbit_
file_map = 'orbit_'

# Image quality of the saved figures.
dpi = 300

# Dictionary orbit_id: number_of_last_orbit for the complete simulation
last_orbits={1: 6000, 2:5000}

# Total number of targets in the sky
total_nb_targets = resx*resy
