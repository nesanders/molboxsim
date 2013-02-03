import matplotlib.image as image
import os,pdb,sys
from pylab import *
from numpy import *
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import proj3d



## Options
#prefix: Prefix for output files
#Temp: Initial temperature
#imageSize: Image size (e.g. '7,5')
#Nf: Nubmer of frames to render
#bps: Ratio of box size to particle size
#N: Number of molecules
#paclph: alpha channel value for rendering particles


if len(sys.argv)>1:  prefix=sys.argv[1]
else: prefix='solid/'

if len(sys.argv)>2:  T=float(sys.argv[2])
else: T=100

if len(sys.argv)>3: imageSize=(float(sys.argv[3].split(',')[0]),float(sys.argv[3].split(',')[1]))
else: imageSize=(8,4.5)

if len(sys.argv)>4:  Nf=int(sys.argv[4])
else: Nf=1000

if len(sys.argv)>5:  bps=float(sys.argv[5])
else: bps=0.5/0.003

if len(sys.argv)>6:  N=int(sys.argv[6])
else: N=200

if len(sys.argv)>7:  palph=float(sys.argv[7])
else: palph=1

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
#N=200 #Number of molecules
bs=5e-1 #box size, meters
ps=bs/bps# 3e-2 # particle size, meters
ts=1e-5 #timestep, s


bounds=[bs,bs,bs]





#Velocities and positions
v=getvel3(T,N,mW)
p=getpos3(bounds,N)
p0=p.copy()


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

os.system('rm '+prefix+'sim*.png')

#Simulate
for n in range(Nf):
  print 'Step',n
  for i in range(3):
    ##Advance one timestep
    p[i]+=v[i]*ts
    ##Use harmonic potential
    v[i]+=-((p[i]-p0[i])/(0.005*ps))**2*sign(p[i]-p0[i])
  ##Add randomness to velocity
  v+=ts*getvel3(T,N,mW)

  #Remake plot
  ax.cla()
  line=ax.scatter(p[0,1:],p[1,1:],p[2,1:],edgecolors='none',s=(imageSize[1]*72*ps/bs)**2,alpha=palph)
  
  # Plot Molly
  x2, y2, z2 = proj3d.proj_transform(p[0,0]*.7,p[1,0]*.7,p[2,0]*.7, ax.get_proj())
  x3,y3=ax.transData.transform((x2, y2))  # convert 2d space to screen space
  im = OffsetImage(im_molly, zoom=0.5)
  ab = AnnotationBbox(im, (x2,y2), xycoords='data', frameon=False)
  ax.add_artist(ab)
  
  fixbounds()
  draw()
  
  plt.show()
  plt.savefig(prefix+'sim'+str(n).zfill(4)+'.png',dpi=160)


os.system('convert -delay 5 -loop 0  -limit memory 256mb -limit map 256mb -quality 100 '+prefix+'sim*.png '+prefix+'animate_sim.mp4')

