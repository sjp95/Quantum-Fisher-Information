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
mpl.rcParams['axes.linewidth'] = 40
mpl.rcParams['axes.edgecolor'] = 'w' 
plt.style.use('classic')

rc('text',usetex=True)
rc('font',family="Times New Roman")

N = 1000#number of points for plotting/interpolation


als=1.0

t=float(sys.argv[1])
t0=int(sys.argv[2])
t1=float(sys.argv[3])
NNN=str(int(sys.argv[4]))
hz=str(int(sys.argv[5]))
ymax=float(sys.argv[6])
uuu=float(sys.argv[7])

plt.figure(figsize=(26,12))
grid = plt.GridSpec(1, 3, wspace=0.5, hspace=0.5,top=0.85, bottom=0.2, right=0.94, left=0.1)
#=============================================================#
                            #DOS#
#=============================================================#
ax = plt.subplot(grid[0, 0])  # Grid for the plot
order=0

s=str(order)
s0=str(t1/10.0)
prename=str("SZ")
if(order==0):
    prename=str("Data/NSC/SZ/SZQW/SZ")
if(order==1):
    prename=str("Data/NSC/SZiSZj/SiSjQW/SZiSiSj")
if(order==2):
    prename=str("Data/NSC/SZiSZjSZk/SZiSjSkQW/SZiSjSk")

ff =str(prename+"_27_100_27_100_"+NNN+"_"+hz+".dat") #"+s+"


file = pb.Path(ff)
#=======================================C-2==========================================#
if file.exists():
    x1 = np.genfromtxt(file, unpack=True)
    x = x1[0, :]
    y = x1[1, :]
    z = abs(x1[t0, :])
    U = z.max() / uuu
    for i in range(len(z)):
        if z[i] > U:
            z[i] = U + np.log10(z[i] / U)

    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')

    c = ax.imshow(zi, extent=[x.min(), x.max(), 0, y.max()], origin='lower', cmap='jet', aspect=5 * (x.max() - x.min()) / y.max(), interpolation='catrom')

    im_ratio = zi.shape[1] / zi.shape[0]

    cb = plt.colorbar(c, ax=ax, fraction=0.045 * im_ratio, pad=0.05, ticks=[min(z), round(max(z) / 2), round(max(z))])
    cb.ax.tick_params(labelsize=45)
    cb.ax.tick_params(width=4, length=7, color='w')
    cb.outline.set_edgecolor('white')
    cb.outline.set_linewidth(4)
    # ax.clim(min(z), round(max(z)))

    plt.ylim(0, ymax)
    plt.xlim(x.min(), x.max())

    plt.xticks([x.min(), (x.max() + x.min()) / 2, x.max()], [r'$0$', r'$\pi$', r'$2\pi$'], fontsize=45)
    plt.yticks([0, 1, 2, 3], fontsize=45)
    ax.tick_params(width=2, length=10, color='w')

    plt.xlabel(r'${\bf q}$', fontsize=45, labelpad=-2)
    plt.ylabel(r'$\frac{\omega}{J_1}$', fontsize=45, labelpad=2)

    
    xx1=(max(x)+min(x))/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    
    yy1=(ymax+min(y))/2
    yy1=(ymax+yy1)/2
    yy1=(ymax+yy1)/2
    #yy1=(ymax+yy1)/2
    plt.text(xx1, yy1, '(a)',fontsize=35,color='white')#-----------------------NO
    
    xx2=(max(x)+min(x))/2
    #xx2=(xx1*2+xx2)
    # xx2=(max(x)+xx2)/2
    # xx2=(max(x)+xx2)/2
    
    plt.text(xx2, yy1, r'$\chi^{NSC}_0({\bf q},\omega)$',fontsize=35,color='white')#-----------------------NO
    
    for spine in ax.spines.values():
        spine.set_edgecolor('w')  # Change frame color
        spine.set_linewidth(4)      # Adjust frame thickness
    #plt.text(xx1, yy1, r'$Order=$ '+s,fontsize=45,color='white')#-----------------------NO
    #plt.text(xx2, yy1, r'$J_x=$ '+s0+r'$J_y$',fontsize=45,color='white')
# #=======================================C-1==========================================#
ax =plt.subplot(grid[0, 1])
order=1

s=str(order)
s0=str(t1/10.0)
prename=str("SZ")
if(order==0):
    prename=str("Data/NSC/SZ/SZQW/SZ")
if(order==1):
    prename=str("Data/NSC/SZiSZj/SiSjQW/SZiSiSj")
if(order==2):
    prename=str("Data/NSC/SZiSZjSZk/SZiSjSkQW/SZiSjSk")

ff =str(prename+"_27_100_27_100_"+NNN+"_"+hz+".dat") #"+s+"


file = pb.Path(ff)
#=======================================C-2==========================================#
if file.exists():
    x1 = np.genfromtxt(file, unpack=True)
    x=x1[0,:]
    y=x1[1,:]
    z=(abs(x1[t0,:]))
    U=z.max()/uuu
    for i in range(len(z)):
        if(z[i]>U):
            z[i]=U+np.log10(z[i]/U)
    print(U)
    #z=z/max(z)
    # z=np.log10(z)
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    c=plt.imshow(zi, extent=[x.min(),x.max(), 0, y.max()], origin='lower',cmap='jet',aspect=5*(x.max()-x.min())/y.max(),interpolation='catrom') #gist_ncar
    im_ratio = zi.shape[1]/zi.shape[0]
    #plt.colorbar(c)
    cb=plt.colorbar(c,fraction=0.045*im_ratio , pad=0.05,ticks=[min(z), round(max(z)/(t))/2,round(max(z)/(t))])#ticks=[0, 0.5,1]
    cb.ax.tick_params(labelsize=45)
    cb.ax.tick_params(width=4,length=7,color='w')
    cb.outline.set_edgecolor('white')
    cb.outline.set_linewidth(4)
    plt.clim(min(z),round(max(z)/(t)))
    #---------------------------------------------------------------------------------------------------------------------------
    #if(order==0):
    xi=np.linspace(x.min(), x.max(),100)
    fx=(np.pi/2.0)*np.sin(abs(xi))
    gx=(np.pi)*np.sin(abs(xi)/2)

    # plt.plot(xi,fx,'--',linewidth=3,color='w')
    # plt.plot(xi,gx,'--',linewidth=3,color='w')
    plt.ylim(0,ymax)
    plt.xlim(x.min(), x.max())

    #---------------------------------------------------------------------------------------------------------------------------
    plt.xticks((x.min(), (x.max()+x.min())/2, x.max()), (r'$0$', r'$\pi$',r'$2\pi$'),fontsize=45)
    plt.yticks([0,1,2,3],fontsize=45)  #[0,2,4,6,8],
    plt.tick_params(width=2,length=10,color='w')
    plt.xlabel(r'${\bf q}$',fontsize=45,labelpad=-2)
    plt.ylabel(r'$\frac{\omega}{J_1}$',fontsize=45,labelpad=2)
    
    xx1=(max(x)+min(x))/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    
    yy1=(ymax+min(y))/2
    yy1=(ymax+yy1)/2
    yy1=(ymax+yy1)/2
    #yy1=(ymax+yy1)/2
    plt.text(xx1, yy1, '(b)',fontsize=35,color='white')#-----------------------NO
    
    xx2=(max(x)+min(x))/2
    #xx2=(xx1*2+xx2)
    # xx2=(max(x)+xx2)/2
    # xx2=(max(x)+xx2)/2
    
    plt.text(xx2, yy1, r'$\chi^{NSC}_1({\bf q},\omega)$',fontsize=35,color='white')#-----------------------NO
    for spine in ax.spines.values():
        spine.set_edgecolor('w')  # Change frame color
        spine.set_linewidth(4)      # Adjust frame thickness
    #plt.text(xx1, yy1, r'$Order=$ '+s,fontsize=45,color='white')#-----------------------NO
    #plt.text(xx2, yy1, r'$J_x=$ '+s0+r'$J_y$',fontsize=45,color='white')
#=======================================C-3==========================================#
plt.subplot(grid[0, 2])
order=2

s=str(order)
s0=str(t1/10.0)
prename=str("SZ")
if(order==0):
    prename=str("Data/NSC/SZ/SZQW/SZ")
if(order==1):
    prename=str("Data/NSC/SZiSZj/SiSjQW/SZiSiSj")
if(order==2):
    prename=str("Data/NSC/SZiSZjSZk/SZiSjSkQW/SZiSjSk")

ff =str(prename+"_27_100_27_100_"+NNN+"_"+hz+".dat") #"+s+"


file = pb.Path(ff)
#=======================================C-2==========================================#
if file.exists():
    x1 = np.genfromtxt(file, unpack=True)
    x=x1[0,:]
    y=x1[1,:]
    z=(abs(x1[t0,:]))
    U=z.max()/uuu
    for i in range(len(z)):
        if(z[i]>U):
            z[i]=U+np.log10(z[i]/U)
    print(U)
    #z=z/max(z)
    # z=np.log10(z)
    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    c=plt.imshow(zi, extent=[x.min(),x.max(), 0, y.max()], origin='lower',cmap='jet',aspect=5*(x.max()-x.min())/y.max(),interpolation='catrom') #gist_ncar
    im_ratio = zi.shape[1]/zi.shape[0]
    #plt.colorbar(c)
    cb=plt.colorbar(c,fraction=0.045*im_ratio , pad=0.05,ticks=[min(z), round(max(z)/(t))/2,round(max(z)/(t))])#ticks=[0, 0.5,1]
    cb.ax.tick_params(labelsize=45)
    cb.ax.tick_params(width=4,length=7,color='w')
    cb.outline.set_edgecolor('white')
    cb.outline.set_linewidth(4)
    plt.clim(min(z),round(max(z)/(t)))
    #---------------------------------------------------------------------------------------------------------------------------
    #if(order==0):
    xi=np.linspace(x.min(), x.max(),100)
    fx=(np.pi/2.0)*np.sin(abs(xi))
    gx=(np.pi)*np.sin(abs(xi)/2)

    # plt.plot(xi,fx,'--',linewidth=3,color='w')
    # plt.plot(xi,gx,'--',linewidth=3,color='w')
    plt.ylim(0,ymax)
    plt.xlim(x.min(), x.max())

    #---------------------------------------------------------------------------------------------------------------------------
    plt.xticks((x.min(), (x.max()+x.min())/2, x.max()), (r'$0$', r'$\pi$',r'$2\pi$'),fontsize=45)
    plt.yticks([0,1,2,3],fontsize=45)  #[0,2,4,6,8],
    plt.tick_params(width=2,length=10,color='w')
    plt.xlabel(r'${\bf q}$',fontsize=45,labelpad=-2)
    plt.ylabel(r'$\frac{\omega}{J_1}$',fontsize=45,labelpad=2)
    
    xx1=(max(x)+min(x))/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    xx1=(min(x)+xx1)/2
    
    yy1=(ymax+min(y))/2
    yy1=(ymax+yy1)/2
    yy1=(ymax+yy1)/2
    #yy1=(ymax+yy1)/2
    plt.text(xx1, yy1, '(c)',fontsize=40,color='white')#-----------------------NO
    
    xx2=(max(x)+min(x))/2
    #xx2=(xx1*2+xx2)
    # xx2=(max(x)+xx2)/2
    # xx2=(max(x)+xx2)/2
    
    plt.text(xx2, yy1, r'$\chi^{NSC}_2({\bf q},\omega)$',fontsize=35,color='white')#-----------------------NO
    #plt.text(xx1, yy1, r'$Order=$ '+s,fontsize=45,color='white')#-----------------------NO
    #plt.text(xx2, yy1, r'$J_x=$ '+s0+r'$J_y$',fontsize=45,color='white')
#=======================================C-1==========================================#





os.makedirs("Plot",exist_ok=True)
ff1 =str("Plot/NSC.pdf")#"+s+"
plt.savefig(ff1,dpi=100)

ff1 =str("Plot/NSC.png")#"+s+"
plt.savefig(ff1,dpi=400)
plt.show()