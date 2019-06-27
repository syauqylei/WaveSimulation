import numpy as np
import matplotlib.pyplot as plt

vel=np.ones([50,50])*4000
np.savetxt('velocity.txt',vel)
for i in range(0,400,10):
    U=np.genfromtxt("output"+str(i)+".txt")
    plt.imshow(U,cmap='gray',vmax=1,vmin=-1)
    plt.show()