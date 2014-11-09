# Resolution of the sky grid
resx = 100
resy = 50


dlambda_over_lambda = 0.8
radius_psf = 15.5
radius_tele = 0.15

# Unit: V magnitude (used for clipping value, important as if mag_max is too high, 
#	the pipeline does not converge)
magnitude_min = 6
magnitude_max = 13

# Unit: ph/ (px s)
sl_min = 1e-8
sl_max = 1e0

# Mirrors efficiency (up to 1 --> 70% = 0.7)
mirror_efficiency_star_class = {'G': 0.65,'K': 0.62,'M': 0.51} 
mirror_efficiency = mirror_efficiency_star_class['G']

# Quantum efficiency for solar type star ie the light we'll get for the stray light
SL_QE=0.67

# Post-treatment stray light reduction (up to 1 --> 99.6% = 0.996)
SL_post_treat_reduction = 0.996

# ppm above the image is considered as contaminated by stray light // USE 5 ppm for bright stars / 70 ppm for faint
ppm_threshold = 10.

# Name of map file without markings. Ie for orbit_000.dat -> orbit_
file_map = 'orbit_'

# Image quality of the saved figures.
dpi = 300

# Dictionary orbit_id: number_of_last_orbit
last_orbits={'700_25_AKTAR': 5322,'800_35_AKTAR_100x50':5210}

# Total number of targets in the sky
total_nb_targets = resx*resy
