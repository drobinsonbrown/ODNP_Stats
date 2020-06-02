import os,sys
import glob
import numpy as np
import matplotlib
#try:
#    matplotlib.use('Agg')
#except ValueError:
#    pass
import matplotlib.pyplot as plt

Ree265 = np.loadtxt(glob.glob('*265*.txt')[0])/10.0
mean265 = np.mean(Ree265)

Ree275 = np.loadtxt(glob.glob('*275*.txt')[0])
mean275 = np.mean(Ree275)

Ree290 = np.loadtxt(glob.glob('*290*.txt')[0])
mean290 = np.mean(Ree290)

Ree310 = np.loadtxt(glob.glob('*310*.txt')[0])
mean310 = np.mean(Ree310)

Ree340 = np.loadtxt(glob.glob('*340*.txt')[0])
mean340 = np.mean(Ree340)

Ree265, bin265 = np.histogram(Ree265,bins=10,density=True)
Ree275, bin275 = np.histogram(Ree275,bins=10,density=True)
Ree290, bin290 = np.histogram(Ree290,bins=10,density=True)
Ree310, bin310 = np.histogram(Ree310,bins=10,density=True)
Ree340, bin340 = np.histogram(Ree340,bins=10,density=True)

l = len(bin265)

colors=['b','r','g','k','c']

plt.plot(bin265[:l-1],Ree265,label='265 K',color=colors[0])
plt.plot(bin275[:l-1],Ree275,label='275 K',color=colors[1])
plt.plot(bin290[:l-1],Ree290,label='290 K',color=colors[2])
plt.plot(bin310[:l-1],Ree310,label='310 K',color=colors[3])
plt.plot(bin340[:l-1],Ree340,label='340 K',color=colors[4])
plt.axvline(x=mean265,color=colors[0],linestyle='--')
plt.axvline(x=mean275,color=colors[1],linestyle='--')
plt.axvline(x=mean290,color=colors[2],linestyle='--')
plt.axvline(x=mean310,color=colors[3],linestyle='--')
plt.axvline(x=mean340,color=colors[4],linestyle='--')

plt.xlabel(r'$R_{ee}$ [nm]')
plt.legend()
plt.savefig('distributions.png')
plt.show()
