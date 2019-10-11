#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import wave.wave as wave
import numpy as np

#create velocty model
vel=np.ones([350,600])*4000
vel[:175,:]=2000

#instantiate wave class
wve = wave.wave(10.0,6,vel)

#Choose the methods "FD2, FD4, FD8, FD12, FD16 or stereo"
wve.set_method('FD16')

#Put a source in a specific location in computation grid.
wve.set_loc(150,300)

#show the result in every n-iteration
wve.set_i_show(250)

#run the simulation(if you choose stereo use run_stereo() instead)
wve.run_fd()

#store the record
rec= wve.record
