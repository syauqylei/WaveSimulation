#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys
fname=str(sys.argv[1])
rec=np.genfromtxt(fname)
plt.pcolormesh(rec,cmap="gray",vmax=0.05,vmin=-0.05)
plt.show()
print(rec.shape)
'''

for i in range(0,6000,100):
    u=np.genfromtxt("output"+str(i)+".txt")
    plt.pcolormesh(u,cmap="gray",vmax=1,vmin=-1)
    plt.show()
'''
