import numpy as np
import matplotlib.pyplot as plt
dt=0.001
nt=100
src=np.zeros([nt,2])
t=np.zeros(nt)
f=50.0

for i in range(nt):
	src[i,1]=(1.-2.*(np.pi*f*(dt*i))**2)*np.exp(-(np.pi*f*(dt*i))**2)
	src[i,0]=dt*i

plt.plot(src[:,0],src[:,1])
plt.show()

np.savetxt("srcfunc.txt",src,fmt='%1.7f')
