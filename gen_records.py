#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import wave.wave as wave
import numpy as np

vel=np.ones([350,600])*4000
vel[:175,:]=2000
wve = wave.wave(10.0,6,vel)

wve.set_method('FD16')
rec=[]

wve.set_loc(150,300)
wve.set_i_show(250)
wve.run_fd()
rec= wve.record