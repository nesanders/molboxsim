

## Maxwell-Boltzmann distribution
kB=1.380648e-23 #J/K
getvel=lambda T,N,M: random.normal(0,sqrt(kB*T/M),N)
getvel3=lambda T,N,M: [getvel(T,N,M),getvel(T,N,M),getvel(T,N,M)]
getpos=lambda bs,N: random.uniform(0,bs,N)
getpos3=lambda bs,N: [getpos(bs,N),getpos(bs,N),getpos(bs,N)]
mW=3e-26 #kg, water molecule


## Initial conditions
N=100 #Number of molecules
T=373 #Temperature, Kelvin
bs=.5 #box size, meters
ts=1e-5 #timestep, s
#Velocities and positions
v=getvel3(T,N,mW)
p=getpos3(bs,N)


## Establish plot
plt.close('all')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
ax.view_init(20,-120)

def fixbounds():
  ax.set_xlim(0,bs)
  ax.set_ylim(0,bs)
  ax.set_zlim(0,bs)
  ax.xaxis.set_ticks([])
  ax.yaxis.set_ticks([])
  ax.zaxis.set_ticks([])

fixbounds()


#Simulate
for n in range(1000):
  for i in range(3):
    p[i]+=v[i]*ts
    ##Which ones hit the edge?
    hit_in=where(p[i]<0)
    hit_out=where(p[i]>bs)
    #Reverse velocity
    v[i][hit_in]=-v[i][hit_in]
    v[i][hit_out]=-v[i][hit_out]
    #Move hits back into box
    p[i][hit_in]=zeros(len(hit_in))
    p[i][hit_out]=ones(len(hit_out))*bs
  ax.cla()
  ax.scatter(p[0],p[1],p[2])
  fixbounds()
  draw()
  #plt.savefig('sim'+str(n)+'.png')


