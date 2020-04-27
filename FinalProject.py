# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 15:25:52 2020

@author: potte
"""

import astropy.io.fits as fits
from astropy.stats import sigma_clip
from astropy import stats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astroML.density_estimation import KNeighborsDensity
from astropy.table import Table

#data from fits file
nsatlas=fits.open('nsa_v0_1_2.fits')
data=nsatlas[1]
mag=22.5-data.data['NMGY']

colors=mag[:,3]-mag[:,4]
redshift=data.data['Z']*10000
#redshift must be made larger; the small numbers (the maximum being 0.055)
#do not work well on the density graph

table=np.column_stack((colors,redshift))
table=np.delete(table,np.argwhere(colors < 0),axis=0)

#Choosing 2000 galaxies at random, making sure they're the same set every time,
#sigma clipping the results, and removing the masked data from the array
N1=2300
np.random.seed(0)
selection=np.random.randint(0,len(table),size=N1)

fit=table[selection]
fit=sigma_clip(fit,sigma=3)
fit=np.delete(fit.data,np.argwhere(fit.mask[:,0] == True),0)

#fit_table=np.delete(fit_table,np.argwhere(fit_table[:,1] < 0),0)
#fit_table=fit_table[(fit_table[:,0] >= 0)]

#Preparing the grid for evaluation. Nx and Ny can be set separately
N=len(fit)
Ny=N
Nx=N
xmin,xmax=(min(fit[:,0]),max(fit[:,0])) 
ymin,ymax=(min(fit[:,1]),max(fit[:,1]))
x=np.linspace(xmin,xmax,Nx)
y=np.linspace(ymin,ymax,Ny)
b=np.vstack(map(np.ravel,np.meshgrid(x,y))).T

#Density estimation fit and evaluation
k=40
knn=KNeighborsDensity('bayesian',k)
knn.fit(fit)
c=knn.eval(b).reshape((Ny,Nx))

#Getting optimization data and science data from what's left
others=np.delete(table.data,selection,0)

N2=1000
np.random.seed(1)
optimize=np.random.randint(0,len(others),size=N2)

opt=others[optimize]
others=np.delete(others,optimize,0)

np.random.seed(2)
science=np.random.randint(0,len(others),size=N2)

sci=others[science]

#Estimating the redshift by finding the highest density corresponding to each
#given color magnitude. Change sci to opt for optimizing
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

sg_idx=np.zeros((N2,1))

for i in range(0,N2):
    sg_idx[i]=find_nearest(x,sci[i,0])

guesses=np.zeros((N2,2))
dens=np.zeros((N2,N))

for i in range(0,N2):
    guesses[i,0]=x[int(sg_idx[i])]
    for j in range(0,N):
        dens[i,j]=c[j,int(sg_idx[i])]
        
densmax=np.zeros((N2,1))        
        
for i in range(0,N2):
    densmax[i]=np.amax(dens[i])
    
h=np.zeros((N2,1))

for i in range(0,N2):
    h,l=np.where(c==densmax[i])
    guesses[i,1]=y[h]

#Results: Plots and difference between accepted and estimated rs
plt.imshow(c,origin='lower',norm=LogNorm(),extent=(xmin,xmax,ymin,ymax),cmap=plt.cm.binary)
plt.title('Density')
plt.xlabel('Color (g-r)')
plt.ylabel('Redshift (z*10000)')
plt.axes().set_aspect('auto')
plt.show()

plt.rcParams.update({'font.size': 15})
plt.scatter(fit[:,0],fit[:,1],s=1,c='k')
plt.scatter(guesses[:,0],guesses[:,1],s=5,c='red')
plt.title('Science Galaxies with Estimated Redshift')
plt.xlabel('Color (g-r)')
plt.ylabel('Redshift (z*10000)')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.show()

#table of differences. Change sci to opt for optimizing
diff=abs(np.around(sci[:,1]-guesses[:,1],decimals=4))
t=Table([np.around(sci[:,1],decimals=4),np.around(guesses[:,1],decimals=4),diff],names=('Accpeted Redshift','Estimated Redshift','Difference'))

d=max(diff)
a=t['Accpeted Redshift'][int(np.where(diff == d)[0])]

print('Maximum Percent Difference: ',(d/a)*100,'%', 'for ',k,' neighbors')