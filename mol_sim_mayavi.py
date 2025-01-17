import matplotlib.image as image
from pylab import *
from numpy import *
from mayavi import mlab
import pdb

## Maxwell-Boltzmann distribution
kB=1.380648e-23 #J/K
getvel=lambda T,N,M: random.normal(0,sqrt(kB*T/M),N)
getvel3=lambda T,N,M: array([getvel(T,N,M),getvel(T,N,M),getvel(T,N,M)])
getpos=lambda bs,N: random.uniform(0,bs,N)
getpos3=lambda bs,N: array([getpos(bs[0],N),getpos(bs[1],N),getpos(bs[2],N)])
mW=3e-26 #kg, water molecule

getdist=lambda p,q: ((p[0]-q[0])**2+(p[1]-q[1])**2+(p[1]-q[1])**2)**.5




## Initial conditions - gas
#N=100 #Number of molecules
#T=373 #Temperature, Kelvin
#bs=.5 #box size, meters
#ps=3e-10 # particle size, meters
#ts=1e-5 #timestep, s

## Initial conditions - liquid
N=100 #Number of molecules
T=20 #Temperature, Kelvin
bs=5e-1 #box size, meters
ps=3e-2 # particle size, meters
ts=1e-5 #timestep, s


## Options
collisions=1
liquid=1
imageSize=(8,4.5)

if liquid: bounds=[bs,bs,bs/2.]
else: bounds=[bs,bs,bs]





#Velocities and positions
v=getvel3(T,N,mW)
p=getpos3(bounds,N)


# read in our png file
im = image.imread('water_mol_s.png')

fig=mlab.figure(bgcolor=(1,1,1),fgcolor=(0.5,0.5,0.5),size=(1280,720))
s=mlab.points3d(p[0],p[1],p[2],colormap='copper',scale_factor=0.05)

#Simulate
@mlab.show
@mlab.animate(delay=500)
def anim():
  ## Establish plot
  fig=mlab.gcf()
  fig.scene.camera.yaw(1)
  #fig.scene.outline()
  for n in range(1000):
    print 'Step',n
    for i in range(3):
      ##Advance one timestep
      p[i]+=v[i]*ts
      ##Which ones hit the edge?
      hit_in=where(p[i]<0)
      hit_out=where(p[i]>bounds[i])
      #Reverse velocity
      v[i][hit_in]=-v[i][hit_in]
      v[i][hit_out]=-v[i][hit_out]
      #Move hits back into box
      p[i][hit_in]=zeros(len(hit_in))
      p[i][hit_out]=ones(len(hit_out))*bounds[i]
    #Search for colliding particles
    if collisions:
      for l in range(N):
        for m in range(N):
          if m!=l and getdist(p[:,m],p[:,l])<ps:
            for i in range(3):
              #Elastic collision equations
              v2f=real(max(roots([2,-2*(v[i,l]+v[i,m]),+2*v[i,l]*v[i,m]])))
              v1f=v[i,l]+v[i,m]-v2f
              v[i,l]=v2f
              v[i,m]=v1f
    ## Plot
    #s.mlab_source.x.setasflat(p[0])
    #s.mlab_source.y.setasflat(p[1])
    #s.mlab_source.z.setasflat(p[2])
    s.mlab_source.set(x=p[0],y=p[1],z=p[2])
    s.scene.save('mav_'+str('n')+'.png')
    yield



anim()
