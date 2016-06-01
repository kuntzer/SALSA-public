import numpy as np
# Resolution of the sky grid
resx = 100
resy = 50

# Aperture of the PSF, in px
aperture_size = 30. # px

# Mirrors efficiency (up to 1 --> 70% = 0.7)
mirror_efficiency_star_class = {'G': 0.65,'K': 0.62,'M': 0.51} 

# Quantum efficiency for solar type star ie the light we'll get for the stray light
SL_QE=0.67

# Post-treatment stray light reduction (up to 1 --> 99.6% = 0.996)
SL_post_treat_reduction = 0.996

# ppm above the image is considered as contaminated by stray light
# USE 5 ppm for bright stars
#    70 ppm for faint
ppm_thresholds = {9: 5, 12: 70}

# Stellar fluxes
# USE G for bright stars
#     K for faint targets
# UNITS: e- / s
stellar_fluxes = {'G':{9: 6.76e5}, # G star flux (Teff=5500) for 330-1100 nm, V=9: total number of photons per second= 1.41x10^6, global throughput = 48%, V=9: total number of electrons per second= 6.76x10^5
		  'K':{12: 4.82e4}, # K star flux (Teff=4500) for 330-1100 nm, V=12: total number of photons per second= 1.1x10^5, global throughput = 44%, V=12: total number of electrons per second= 4.82x10^4
		}

# Name of map file without markings. Ie for orbit_000.dat -> orbit_
file_map = 'orbit_'

# Image quality of the saved figures.
dpi = 300

# Dictionary orbit_id: number_of_last_orbit
last_orbits={'650_5_conf4e':5375, '700_5_conf4e': 5316}

# Total number of targets in the sky
total_nb_targets = resx*resy
