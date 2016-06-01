import numpy as np
import pylab as plt
import os

import resources.figures as figures

###############################################################################

# UNCOMMENT FOR 80%
# A list of dictionaries containing all data to plot.
# Mandatory keys:
# 'folder', 'fname'
# Optional keys:
# "label",'lw', 'color'
"""skycoverage=[
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_0.0pst_0.996SLreduction_cumul_.txt",
	# 'color':"b","label":"No SAA",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	# 'color':"k","label":r"$25^\circ$, no SL",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.996SLreduction_cumul_.txt",
	# 'color':"r","label":"$25^\circ,\ 99.6\%$",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"400 nm",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.980SLreduction_cumul_.txt",
	# 'color':"grey","label":"$25^\circ,\ 98\%$",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.970SLreduction_cumul_.txt",
	 #'color':"maroon","label":"$25^\circ,\ 97\%$",'lw':2},	
	 {'folder':"700_25_conf3_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"conf3",'lw':2},
	 {'folder':"6pm_700_25_400nm_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"6pm",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_79min_V12.0_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"6am / No SAA",'lw':2},
	 {'folder':"6pm_700_25_400nm_misc",
	 'fname':"skycoverage_79min_V12.0_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"6pm / No SAA",'lw':2},	
	]"""

skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 pm",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 pm",'lw':2},
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_80min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, 6 am",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_80min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, 6 pm",'lw':2},	
	]

skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 pm",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 pm",'lw':2},
	]

skycoverage=[
	{'folder':"6pm_650_25_conf3_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, conf3",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, conf4",'lw':2},
	{'folder':"6pm_700_25_conf3_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, conf3",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, conf4",'lw':2},
	{'folder':"6pm_800_25_conf3_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, conf3",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, conf4",'lw':2},	
	]

skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 am",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 am",'lw':2},
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, 6 am",'lw':2},	
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 pm",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 pm",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, 6 pm",'lw':2},	
	]

skycoverage=[
	{'folder':"6am_650_5_conf4e_misc",
	 'fname':"skycoverage_78min_V12.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6am_700_5_conf4e_misc",
	 'fname':"skycoverage_79min_V12.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},	
	]

# What is the target ? If None, nothing is shown
requirement_days=13
requirement_obs=25

# Force the axis, if set to None, automatic from pylab
xlim=[0,60]
ylim=[0,45]

# Filename for the figure, if None, the figure is not saved. Without extension!
out_fname="skycoverage_stat_scireq2.2"
'''
# A list of dictionaries containing all data to plot.
# Mandatory keys:
# 'folder', 'fname'
# Optional keys:
# "label",'lw', 'color'
"""
skycoverage=[
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"b","label":"No SAA",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"k","label":r"$25^\circ$, no SL",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"r","label":"$25^\circ,\ 99.6\%$",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.950SLreduction_cumul_.txt",
	 'color':"g","label":"$25^\circ,\ 95\%$",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.900SLreduction_cumul_.txt",
	 'color':"indigo","label":"$25^\circ,\ 90\%$",'lw':2},
	]

skycoverage=[
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_81min_V12.0_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"b","label":"No SAA",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"k","label":r"$25^\circ$, no SL",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"r","label":"$25^\circ,\ 99.6\%$",'lw':2},
	#{'folder':"800_25_400nm_misc",
	# 'fname':"skycoverage_81min_V12.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"gold","label":"$25^\circ,\ 99\%$",'lw':2},
	#{'folder':"800_25_400nm_misc",
	# 'fname':"skycoverage_81min_V12.0_SAA_0.980SLreduction_cumul_.txt",
	# 'color':"grey","label":"$25^\circ,\ 98\%$",'lw':2},
	#{'folder':"800_25_400nm_misc",
	# 'fname':"skycoverage_81min_V12.0_SAA_0.970SLreduction_cumul_.txt",
	# 'color':"maroon","label":"$25^\circ,\ 97\%$",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.960SLreduction_cumul_.txt",
	 'color':"Navy","label":"$25^\circ,\ 96\%$",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_81min_V12.0_SAA_0.950SLreduction_cumul_.txt",
	 'color':"g","label":"$25^\circ,\ 95\%$",'lw':2},
	]

skycoverage=[
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_50min_V9.0_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"b","label":"No SAA",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	 'color':"k","label":r"$25^\circ$, no SL",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"r","label":"$25^\circ,\ 99.6\%$",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.950SLreduction_cumul_.txt",
	 'color':"g","label":"$25^\circ,\ 95\%$",'lw':2},
	{'folder':"800_25_400nm_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.900SLreduction_cumul_.txt",
	 'color':"indigo","label":"$25^\circ,\ 90\%$",'lw':2},
	]
	
skycoverage=[
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_0.0pst_0.996SLreduction_cumul_.txt",
	# 'color':"b","label":"No SAA",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.0pst_0.996SLreduction_cumul_.txt",
	# 'color':"k","label":r"$25^\circ$, no SL",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.996SLreduction_cumul_.txt",
	# 'color':"r","label":"$25^\circ,\ 99.6\%$",'lw':2},
	{'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"400 nm",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.980SLreduction_cumul_.txt",
	# 'color':"grey","label":"$25^\circ,\ 98\%$",'lw':2},
	#{'folder':"700_25_400nm_misc",
	# 'fname':"skycoverage_79min_V12.0_SAA_0.970SLreduction_cumul_.txt",
	 #'color':"maroon","label":"$25^\circ,\ 97\%$",'lw':2},	
	{'folder':'700_25_conf3_misc',
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"conf3",'lw':2},
	 {'folder':"6pm_700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"6pm",'lw':2},
	 {'folder':"700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"6am / No SAA",'lw':2},
	 {'folder':"6pm_700_25_400nm_misc",
	 'fname':"skycoverage_49min_V9.0_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"6pm / No SAA",'lw':2},

	]
skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 pm",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 pm",'lw':2},
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, 6 am",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, 6 pm",'lw':2},	
	]
skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 pm",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 pm",'lw':2},
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, 6 am",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, 6 pm",'lw':2},	
	]

skycoverage=[
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 pm",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_0.0pst_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 pm",'lw':2},
	]

skycoverage=[
	{'folder':"6am_650_25_conf3_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, conf3",'lw':2},
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, conf4",'lw':2},
	{'folder':"6am_700_25_conf3_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, conf3",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, conf4",'lw':2},
	{'folder':"6am_800_25_conf3_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, conf3",'lw':2},	
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, conf4",'lw':2},	
	]
"""
skycoverage=[
	#{'folder':"6am_650_25_conf3_misc",
	# 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"gold","label":"650 km, conf3",'lw':2},
	{'folder':"6am_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"650 km, 6 am",'lw':2},
	#{'folder':"6am_700_25_conf3_misc",
	# 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"r","label":"700 km, conf3",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"green","label":"700 km, 6 am",'lw':2},
	#{'folder':"6am_800_25_conf3_misc",
	# 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"k","label":"800 km, conf3",'lw':2},	
	{'folder':"6am_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"Navy","label":"800 km, 6 am",'lw':2},
	#{'folder':"6pm_650_25_conf3_misc",
	# 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"gold","label":"650 km, conf3",'lw':2},
	{'folder':"6pm_650_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 pm",'lw':2},
	#{'folder':"6pm_700_25_conf3_misc",
	# 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"r","label":"700 km, conf3",'lw':2},
	{'folder':"6pm_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 pm",'lw':2},
	#{'folder':"6pm_800_25_conf3_misc",
	# 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	# 'color':"k","label":"800 km, conf3",'lw':2},	
	{'folder':"6pm_800_25_conf4_misc",
	 'fname':"skycoverage_50min_V9.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"k","label":"800 km, 6 pm",'lw':2},	
	]
"""
skycoverage=[
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_79min_V14.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"r","label":"80\%",'lw':2},
	{'folder':"6am_700_25_conf4_misc",
	 'fname':"skycoverage_49min_V14.0_SAA_0.990SLreduction_cumul_.txt",
	 'color':"b","label":"50\%",'lw':2},	
	]"""

skycoverage=[
	{'folder':"6am_650_5_conf4e_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"gold","label":"650 km, 6 am",'lw':2},
	{'folder':"6am_700_5_conf4e_misc",
	 'fname':"skycoverage_49min_V9.0_SAA_0.996SLreduction_cumul_.txt",
	 'color':"r","label":"700 km, 6 am",'lw':2},	
	]
# What is the target ? If None, nothing is shown
#requirement_days=13
#requirement_obs=25

requirement_days=50
requirement_obs=50

# Force the axis, if set to None, automatic from pylab
xlim=[0,90]
ylim=[0,90]

# Filename for the figure, if None, the figure is not saved. Without extension!
#out_fname="800_25_400nm_figures/skycoverage_stat_81min"
out_fname="skycoverage_stat_scireq2.1"
'''
# Nicer plot + exportation in pdf, eps
fancy=True

###############################################################################

if fancy:figures.set_fancy()
legend=False
f=plt.figure()
for el in skycoverage:
	data=np.loadtxt(os.path.join(el['folder'],el['fname']))
	if 'lw' in el:
		lw=el['lw']
	else:
		lw=2
	if 'label' in el:
		label=el['label']
		legend=True
	else:
		label=None
	if 'color' in el:
		color=el['color']

		plt.plot(data[:,0],data[:,1],lw=lw,label=label,color=color)
	else:
		plt.plot(data[:,0],data[:,1],lw=lw,label=label)


if legend:
	plt.legend(loc='best')

if not requirement_days is None:
	plt.axvline(x=requirement_days, c="gray", lw=2)
if not requirement_obs is None:
	plt.axhline(y=requirement_obs, c="gray", lw=2)

plt.ylabel(r"$\mathrm{Fraction\ of\ full\ sky\ [\%]}$")
plt.xlabel(r"$\mathrm{Required\ observation\ time\ [days]}$")

if not xlim is None:
	plt.xlim(xlim)
if not ylim is None:
	plt.ylim(ylim)
plt.grid()

if not out_fname is None:
	figures.savefig(out_fname, f, fancy)
	print 'saved as %s' % out_fname
plt.show()
