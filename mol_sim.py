import matplotlib.image as image
from pylab import *
from numpy import *
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import proj3d

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
N=200 #Number of molecules
T=100 #Temperature, Kelvin
bs=5e-1 #box size, meters
ps=3e-2 # particle size, meters
ts=1e-5 #timestep, s


## Options
collisions=1
liquid=1
prefix='liquid/'
imageSize=(8,4.5)

if liquid: bounds=[bs,bs,bs/2.]
else: bounds=[bs,bs,bs]





#Velocities and positions
v=getvel3(T,N,mW)
p=getpos3(bounds,N)


## Establish plot
plt.close('all')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=imageSize)
ax = Axes3D(fig)
ax.view_init(20,-120)


def draw_l_surf():
  """Plot water surface"""
  verts=[list(zip([0,bs,bs,0],[0,0,bs,bs]))]
  poly = PolyCollection(verts, facecolors = 'b')
  poly.set_alpha(0.2)
  ax.add_collection3d(poly, zs=[bs/2.+ps/2.], zdir='z')


def fixbounds():
  ax.set_xlim(0,bs)
  ax.set_ylim(0,bs)
  ax.set_zlim(0,bs)
  ax.xaxis.set_ticks([])
  ax.yaxis.set_ticks([])
  ax.zaxis.set_ticks([])

fixbounds()
## Get data bounds
xmin, ymin, zmin = proj3d.proj_transform(bs,bs,bs/2.+ps/2., ax.get_proj())
xmax, ymax, zmax = proj3d.proj_transform(0,0,0, ax.get_proj())

# read in our png file
im_mod = image.imread('water_mol_s.png')
im_molly = image.imread('Molly1_sm.png')


#Simulate
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
      for m in range(1,N/2):
        if getdist(p[:,m],p[:,l])<ps:
          for i in range(3):
            #Elastic collision equations
            v2f=real(max(roots([2,-2*(v[i,l]+v[i,m]),+2*v[i,l]*v[i,m]])))
            v1f=v[i,l]+v[i,m]-v2f
            v[i,l]=v2f
            v[i,m]=v1f
  ax.cla()
  line=ax.scatter(p[0,1:],p[1,1:],p[2,1:],edgecolors='none',s=(imageSize[1]*72*ps/bs)**2)
  
  if liquid: draw_l_surf()
  
  # Plot Molly
  x2, y2, z2 = proj3d.proj_transform(p[0,0]-ps*5,p[1,0]-ps*5,p[2,0]-ps*5, ax.get_proj())
  x3,y3=ax.transData.transform((x2, y2))  # convert 2d space to screen space
  im = OffsetImage(im_molly, zoom=0.5)
  ab = AnnotationBbox(im, (x2,y2), xycoords='data', frameon=False)
  ax.add_artist(ab)
  
  fixbounds()
  draw()
  
  plt.show()
  plt.savefig(prefix+'sim'+str(n).zfill(4)+'.png',dpi=160)


