import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import pathlib as pb
import os
import glob
import sys
from matplotlib import rc
import matplotlib.pyplot as mpl
mpl.rcParams['font.size'] = 30
mpl.rcParams['axes.linewidth'] = 10
plt.style.use('classic')

rc('text',usetex=True)
rc('font',family="Times New Roman")
rc('axes', linewidth=5)

N = 1000#number of points for plotting/interpolation

chop=10
als=1.0
plt.figure(figsize=(20,20))
grid = plt.GridSpec(1, 1, wspace=0.05, hspace=0.02, bottom=0.2)
#=============================================================#
                            #DOS#
#=============================================================#

#=======================================C-1==========================================#
plt.subplot(grid[0, 0])
t=float(sys.argv[1])
t0=int(sys.argv[2])
#
t1=float(sys.argv[3])
order=int(sys.argv[4])
s=str(order)
s0=str(t1/10.0)

prename=str("SZ")
if(order==0):
    prename=str("Data/NSC/SZ/SZQW/SZ")
if(order==1):
    prename=str("Data/SC/SiSj/SiSjQW/SiSj")
if(order==2):
    prename=str("SjSk")

ff =str(prename+"_0_0_100_100_36_18_18.dat") #"+s+"

file = pb.Path(ff)
#=======================================C-2==========================================#
if file.exists():
    x1 = np.genfromtxt(file, unpack=True)
    x=x1[0,:]
    y=x1[1,:]
    z=abs(x1[t0,:])
    z=z/max(z)
    # z=np.log10(-x1[3,:]+1)
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='nearest')
    c=plt.imshow(zi, extent=[x.min(),x.max(), 0, y.max()], origin='lower',cmap='magma',aspect=(x.max()-x.min())/y.max(),interpolation='catrom')
    im_ratio = zi.shape[1]/zi.shape[0]
    #plt.colorbar(c)
    cb=plt.colorbar(c,fraction=0.045*im_ratio , pad=0.05,ticks=[min(z), max(z)/(2*t),max(z)/(t)])#ticks=[0, 0.5,1]
    cb.ax.tick_params(labelsize=60)
    cb.ax.tick_params(width=4,length=7)
    plt.clim(min(z),max(z)/(t))
    #---------------------------------------------------------------------------------------------------------------------------
    if(order==0):
        xi=np.linspace(x.min(), x.max(),100)
        fx=(np.pi/2.0)*np.sin(abs(xi))
        gx=(np.pi)*np.sin(abs(xi)/2)

        plt.plot(xi,fx,'--',linewidth=3,color='w')
        plt.plot(xi,gx,'--',linewidth=3,color='w')
        plt.ylim(0,y.max())
        plt.xlim(x.min(), x.max())

    #---------------------------------------------------------------------------------------------------------------------------
    plt.xticks((x.min(), 0, x.max()), (r'$-\pi$', r'$0$',r'$\pi$'),fontsize=45)
    plt.yticks(fontsize=45) #[0,2,4,6,8],
    plt.tick_params(width=4,length=10)
    plt.xlabel(r'${\bf q}$',fontsize=45,labelpad=-2)
    plt.ylabel(r'$\frac{\omega}{J_1}$',fontsize=45,labelpad=2)
    plt.text(min(x)+24*(max(x)-min(x))/36, max(y)-1*max(y)/8, r'$Order=$ '+s,fontsize=45,color='white')#-----------------------NO
    plt.text(min(x)+12*(max(x)-min(x))/36, max(y)-1*max(y)/8, r'$J_2=$ '+s0+r'$J_1$',fontsize=45,color='white')
#=======================================C-1==========================================#
#=======================================C-2==========================================#
# N = 100
# plt.subplot(grid[0, 1])
# t=int(10)
# s=str(t)

# t1=float(t)
# s0=str(t1/100.0)

# ff =str("Data/Twoc_1.dat") #"+s+"
# file = pb.Path(ff)
# #=======================================C-2==========================================#
# if file.exists():
#     x1 = np.genfromtxt(file, unpack=True)
#     x=x1[0,:]
#     y=x1[1,:]
#     z=x1[2,:]
#     xi = np.linspace(x.min(), x.max(), N)
#     yi = np.linspace(y.min(), y.max(), N)
#     zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
#     c=plt.imshow(zi, extent=[0,x.max(), 0, y.max()], origin='lower',cmap='magma',aspect=x.max()/y.max(),interpolation='catrom')
#     im_ratio = zi.shape[1]/zi.shape[0]
#     # cb=plt.colorbar(c,fraction=0.045*im_ratio , pad=0.05,ticks=[0, 0.5,1])
#     # cb.ax.tick_params(labelsize=60)
#     # cb.ax.tick_params(width=4,length=7)
#     plt.clim(0,1)
#     plt.xticks((0, 2, 4,6), ('$(0,0)$', '$(\pi, 0)$','$(\pi, \pi)$', '$(0,0)$'),fontsize=45)
#     plt.yticks([0,2,4,6,8],fontsize=45)
#     plt.tick_params(width=4,length=10)
#     plt.xlabel(r'${\bf q}$',fontsize=45,labelpad=-2)
#     plt.ylabel(r'$\frac{\omega}{J_1}$',fontsize=45,labelpad=2)
#     plt.text(max(x)-1.5*max(x)/8, max(y)-1*max(y)/8, '(b)',fontsize=45,color='white')#-----------------------NO
#     plt.text(max(x)/8, max(y)-1*max(y)/8, r'$J_2=$ '+s0+r'$J_1$',fontsize=45,color='white')
# #=======================================C-2==========================================#
# #=======================================C-3==========================================#
# plt.subplot(grid[0, 2])
# t=int(10)
# s=str(t)

# t1=float(t)
# s0=str(t1/100.0)

# ff =str("Data/Twoc_2.dat") #"+s+"
# file = pb.Path(ff)
# #=======================================C-3==========================================#
# if file.exists():
#     x1 = np.genfromtxt(file, unpack=True)
#     x=x1[0,:]
#     y=x1[1,:]
#     z=x1[2,:]
#     xi = np.linspace(x.min(), x.max(), N)
#     yi = np.linspace(y.min(), y.max(), N)
#     zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
#     c=plt.imshow(zi, extent=[0,x.max(), 0, y.max()], origin='lower',cmap='magma',aspect=x.max()/y.max(),interpolation='catrom')
#     im_ratio = zi.shape[1]/zi.shape[0]
#     # cb=plt.colorbar(c,fraction=0.045*im_ratio , pad=0.05)#,ticks=[0, 0.5,1]
#     # cb.set_ticks([0, 1])
#     # cb.set_ticklabels(['Min', 'Max'])
#     # cb.ax.tick_params(labelsize=60)
#     # cb.ax.tick_params(width=4,length=7)
#     plt.clim(0,1)
#     plt.xticks((0, 2, 4,6), ('$(0,0)$', '$(\pi, 0)$','$(\pi, \pi)$', '$(0,0)$'),fontsize=45)
#     plt.yticks([0,2,4,6,8],fontsize=45)
#     plt.tick_params(width=4,length=10)
#     plt.xlabel(r'${\bf q}$',fontsize=45,labelpad=-2)
#     plt.ylabel(r'$\frac{\omega}{J_1}$',fontsize=45,labelpad=2)
#     plt.text(max(x)-1.5*max(x)/8, max(y)-1*max(y)/8, '(c)',fontsize=45,color='white')#-----------------------NO
#     plt.text(max(x)/8, max(y)-1*max(y)/8, r'$J_2=$ '+s0+r'$J_1$',fontsize=45,color='white')
# #=======================================C-3==========================================#




os.makedirs("Plot/SC",exist_ok=True)
ff1 =str("Plot/SC/Order_"+s+".pdf")#
plt.savefig(ff1,dpi=100)
plt.show()