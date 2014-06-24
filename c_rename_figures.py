''' copies files to tmp_figures/ to prepare making the movie (files must have a continous naming, starting from %07d % 0 '''

import os

repeat = 2

# March 21
t_ini = 114000
t_end = 118023

t_ini = 1
t_end = 26

# test
t_ini = 247848
t_end = 247948


t_ini = 290215 # orbit 2878
t_end = 290315 # orbit 2878

t_ini = 202959 # orbit 2013
t_end = 203059 # orbit 2013

t_ini = 247142 # orbit 2451
t_end = 247242 # orbit 2451

t_ini = 510321 # orbit 5060
t_end = 510421 # orbit 5060
#t_ini = 4036
#t_end = 4136 

minute = t_ini
iterator = 0
r = 1
while (minute < t_end):
#	fname_in = 'orbitID_figures/debug/dist_earthlimb_%07d.png' % (minute)
#	fname_out = 'figures_tmp/dist_earthlimb_%07d.png' % (iterator)

	fname_in = 'orbitID_figures/debug/flux_%07d.png' % (minute)
	fname_out = 'figures_tmp/flux_%07d.png' % (iterator)

#	fname_in = 'orbitID_flux/flux_%d.dat' % (minute)
#	fname_out = 'figures_tmp/flux_%d.dat' % (iterator)

	command = "cp %s %s" % (fname_in, fname_out)
	os.system(command)
	if r == repeat:
		minute += 1
		r = 1
	else: r+=1
	iterator += 1
