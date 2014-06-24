######################################################################
# DEFINITIONS AND INCLUDES
import numpy as np
import subprocess

from resources.routines import *
import parameters as param
######################################################################
# PARAMETERS

movie_name = 'straylight_orbitID_o1.mkv'

subprocess.call(["avconv","-f","image2","-i","figures_tmp/flux_%07d.png","-c:v","libx264","-r","240","-s","hd720","-crf","16",movie_name])

#subprocess.call(["avconv","-f","image2","-i","figures_tmp/flux_%07d.png","-c:v","libx264","-s","hd720","-crf","3",movie_name])
