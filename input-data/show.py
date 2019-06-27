import numpy as np
import matplotlib.pyplot as plt
#vel=np.ones([100,100])*4000.0
#np.savetxt("velocity.txt",vel)
for i in range(0,2000,10):
    U=np.genfromtxt("output"+str(i)+".txt")
    plt.imshow(U,cmap='gray',vmax=0.1,vmin=-0.1)
    plt.show()

U=np.genfromtxt("rec.txt")
plt.imshow(U,cmap='gray',vmax=0.1,vmin=-0.1)
plt.show()
