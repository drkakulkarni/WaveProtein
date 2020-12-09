#!/usr/bin/env python
import matplotlib.pyplot as plt
import csv
import os
import numpy as np
import pandas as pd
import scipy as sp
import sys
import pycwt as wavelet
from pycwt.helpers import find
from matplotlib.image import NonUniformImage
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
from numpy import matrix
from tabulate import tabulate


###################################################################################################################################
#PLEASE NOTE: YOU HAVE TO EDIT THE PATH OF WCT SIGNIFICANCE FILE PATH IN " site-packages/pycwt/wavelets.py" (Line numbers 568 & 632).
#IT SHOULD BE SAME AS IN THE CURRENT FILE, LINE NUMEBRS 132, 135 AND 215
################################################################################################################################## 

while True:
	sys.stdout.write('Have you edited this script for Prptein_name start and end residue numbers of P-loop, SWI, SWII regions?\n')
	query = raw_input().lower()
	if query == "yes":
		dum =1
		#print('Please answer with yes or no!')
		break
	if query == 'no':
		print("please edit the start and end residue numbers of P-loop, SWI, SWII regions\n") 
		exit()
	else:
		print(" In valid answer. Plase input yes or no\n")

###################################################################################################
Protein_name ="Cdc42"
SWI_start = 25
SWI_end = 39

SWII_start = 56
SWII_end = 66

Ploop_start =9
Ploop_end = 16
###################################################################################################

SMALL_SIZE = 18
MEDIUM_SIZE = 19
BIGGER_SIZE = 22


font_kk = {'family': 'arial',
        'color':  'black',
        'weight': 'bold',
        'size': 20,
        }

plt.rc('font', size=SMALL_SIZE,weight='bold')          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE )    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



#scg.set_default_wavelet('cmor1-1.5')
#scg.set_default_wavelet('morl')

#############################################Read Files###################################
root_name1 = sys.argv[1]

file_name1=root_name1+".csv"


if os.path.exists(file_name1):
	dum=1
else:
  print("The file %s does not exist" % file_name1)
  exit()

col_names=['Res', 'Rnum', 'Rco']
df = pd.read_csv(file_name1, names=col_names)

R1=df['Res'].astype(str)
df['Rnum']=df['Rnum'].astype('int64')
df['Rco']=df['Rco'].astype('float64')

t1 = df['Rnum'].to_numpy() 
s1 = df['Rco'].to_numpy() 

#title = 'Wavelet'
label = 'RCV '
units = 'Res No'

root_name2 = sys.argv[2]

file_name2=root_name2+".csv"

if os.path.exists(file_name2):
	dum=1
else:
  print("The file %s does not exist" % file_name2)
  exit()

df = pd.read_csv(file_name2, names=col_names)

R2=df['Res'].astype(str)
df['Rnum']=df['Rnum'].astype('int64')
df['Rco']=df['Rco'].astype('float64')

t2 = df['Rnum'].to_numpy() 
s2 = df['Rco'].to_numpy() 

print("\n")
if (np.array_equal(t1,t2)):
	print("The two signals are of identical lengths")
else:
	comb_name = file_name1+","+file_name2
	print("CSV files have diffrent lenths %s " % comb_name)
	exit()

for i, val in enumerate(R1):
	if (R1[i] != R2[i]):
		print("Res no {0} has diffrenet residues {1} in {2} and {3} in {4}" .format( t1[i], R1[i],root_name1, R2[i], root_name2))
##################################################################################################################

n1 = t1.size
n2 = t2.size
n = min(n1, n2)

s1 = s1/s1.std()
s2 = s2/s2.std()


file_rot = '/Users/admin/.cache/'
sim_file = file_rot+root_name1+"_"+root_name2+".txt"
if os.path.exists(sim_file):
  os.rename(sim_file , '/Users/admin/.cache/pycwt/kl.txt')
else:
  print("The WC significance file does not exist")

#########################################################Wavelet Transform of signal 1###########

mother = wavelet.Morlet(6)          # Morlet mother wavelet with m=6
slevel = 0.95                       # Significance level
dj = 1.0/12.0                           # Twelve sub-octaves per octaves
s0 = -1  # 2 * dt                   # Starting scale
J = -1  # 7 / dj                    # Seven powers of two with dj sub-octaves
dt = 1
if True:
    alpha1, _, _ = wavelet.ar1(s1)  # Lag-1 autocorrelation for red noise
    alpha2, _, _ = wavelet.ar1(s2)  # Lag-1 autocorrelation for red noise
else:
    alpha1 = alpha2 = 0.0           # Lag-1 autocorrelation for white noise

# The following routines perform the wavelet transform and siginificance
# analysis for two data sets.
W1, scales1, freqs1, coi1, _, _ = wavelet.cwt(s1, dt, dj, s0, J, mother)
signif1, fft_theor1 = wavelet.significance(1.0, dt, scales1, 0, alpha1, significance_level=0.95, wavelet=mother)

####### Calculations#########################################################
power1 = (np.abs(W1)) ** 2             # Normalized wavelet power spectrum
period1 = 1/freqs1
sig95_1 = np.ones([1, n1]) * signif1[:, None]
sig95_1 = power1 / sig95_1             # Where ratio > 1, power is significant

#################################################################################################

#########################################################Wavelet Transform of signal 2###########

W2, scales2, freqs2, coi2, _, _ = wavelet.cwt(s2, dt, dj, s0, J, mother)
signif2, fft_theor2 = wavelet.significance(1.0, dt, scales2, 0, alpha2, significance_level=0.95, wavelet=mother)

####### Calculations#########################################################
power2 = (np.abs(W2)) ** 2             # Normalized wavelet power spectrum
period2 = 1/freqs2
sig95_2 = np.ones([1, n2]) * signif1[:, None]
sig95_2 = power2 / sig95_2             # Where ratio > 1, power is significant

#################################################################################################

######################################################### Cross Wavelet Transform of signal 1 &  2 ###########

#s2, _, _ = wavelet.helpers.boxpdf(s2)
#s1, _, _ = wavelet.helpers.boxpdf(s1)
#s2 = (s2-s2.mean())/s2.std()
#s1 = (s1-s1.mean())/s1.std()

XWT,  coi_cr, freqs_cr, signif_cr = wavelet.xwt(s1, s2, 1, dj=dj, s0=s0, J=-1, significance_level=0.95, wavelet=mother, normalize=True)

####### Calculations#########################################################
power_cr = (np.abs(XWT)) ** 2             # Normalized wavelet power spectrum
sig_cr = np.ones([1, n]) * signif_cr[:, None]
period_cr = 1/freqs_cr
sig95_cr = power_cr / sig_cr             # Where ratio > 1, power is significant

#################################################################################################

######################################################### Wavelet Coharance signal 1 &  2 ###########

WCT, aWCT, cor_coi, freq, sig = wavelet.wct(s1, s2, dt, dj=dj, s0=-1, J=-1,
                                             significance_level=0.8646,
                                             wavelet=mother, normalize=True,
                                             cache=True)
WCT[WCT > 1.0] = 1.0
WCT[WCT < -1.0] = -1.0
####### Calculations#########################################################
cor_power = np.abs(WCT)
cor_sig = np.ones([1, n]) * sig[:, None]
cor_sig = np.abs(WCT) / cor_sig  # Power is significant where ratio > 1
cor_period = 1 / freq

angle = 0.5 * np.pi - aWCT
u, v = np.cos(angle), np.sin(angle)

KK_cutoff=float(sys.argv[3])
T_idx = (np.abs(np.log2(cor_period)-KK_cutoff)).argmin()
print("Long 2 T is {0} and Frequeny is {1} and index is {2}".format(np.log2(cor_period[T_idx]),cor_period[T_idx], T_idx))

print(np.shape(u))

#################################################################################################

os.rename('/Users/admin/.cache/pycwt/kl.txt', sim_file )

print("P-Loop Start= {0} {1} and End= {2} {3}" .format(t1[Ploop_start-t1[0]], R1[Ploop_start-t1[0]], t1[Ploop_end-t1[0]], R1[Ploop_end-t1[0]]))
print("SWI Start= {0} {1} and End= {2} {3}" .format(t1[SWI_start-t1[0]], R1[SWI_start-t1[0]], t1[SWI_end-t1[0]], R1[SWI_end-t1[0]]))
print("SWII Start= {0} {1} and End= {2} {3}\n" .format(t1[SWII_start-t1[0]], R1[SWII_start-t1[0]], t1[SWII_end-t1[0]], R1[SWII_end-t1[0]]))


########################Calculate the angle between vectors########################
#def get_anglesum(u,v,nxy,a1,a2):
#	b=int(sys.argv[3])
#	U1=np.array(u[:b:1, a1:a2:1])
#	U2=np.sum(U1, axis=1)
#	U3=np.sum(U2)
#
#	V1=np.array(v[:b:1, a1:a2:1])
#	V2=np.sum(V1, axis=1)
#	V3=np.sum(V2)
#	Angle=np.degrees(np.arctan(V3/U3))
#	if(Angle<1):
#		Angle=Angle*(-1)
#	return(Angle, U3/(b*(a2-a1)), V3/(b*(a2-a1)))

def get_wanglesum(u,v,nxy,wy,a1,a2):
	b=int(T_idx)+nxy
	U1=np.array(np.multiply(u[:b, a1:a2] , wy[:b, a1:a2]))
	U2=np.sum(U1, axis=1)
	U3=np.sum(U2)

	V1=np.array(np.multiply(v[:b, a1:a2] , wy[:b, a1:a2]))
	V2=np.sum(V1, axis=1)
	V3=np.sum(V2)
	Angle=np.degrees(np.arctan(V3/U3))
	if(Angle<1):
		Angle=Angle*(-1)
	#return(Angle, U3/(b*(a2-a1)), V3/(b*(a2-a1)))
	return(Angle, U3, V3)
#np.sin(np.deg2rad(90))

P_loop_ang, pu,pv=get_wanglesum(u,v,0,cor_power,Ploop_start-t1[0], Ploop_end-t1[0])
SWI_ang,s1u,s1v=get_wanglesum(u,v,20,cor_power,SWI_start-t1[0], SWI_end-t1[0])
SWII_ang,s2u,s2v=get_wanglesum(u,v,20,cor_power,SWII_start-t1[0], SWII_end-t1[0])

#P_loop_ang, pu,pv=get_anglesum(u,v,0,Ploop_start-t1[0], Ploop_end-t1[0])
#SWI_ang,s1u,s1v=get_anglesum(u,v,20,SWI_start-t1[0], SWI_end-t1[0])
#SWII_ang,s2u,s2v=get_anglesum(u,v,30,SWII_start-t1[0], SWII_end-t1[0])

print("P-ang={0}, pu={1}, pv={2}".format(P_loop_ang, pu,pv))
print("SW1-ang={0}, s1u={1}, s1v={2}".format(SWI_ang, s1u,s1v))
print("P-ang={0}, s2u={1}, s2v={2}\n".format(SWII_ang, s2u,s2v))


########################Get angle between vector and Y-axis###########################

def Y_ang(ak, bk,u,v):
	vector_1 = np.array([ak,bk])
	vector_2 = np.array([u, v])

	unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
	unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
	dot_product = np.dot(unit_vector_1, unit_vector_2)
	AG= np.arccos(dot_product)
	return(np.degrees(AG))


######################################Plot the Signal(RCO)##############################
plt.ioff()
figprops = dict(num=1, figsize=(11, 8), dpi=100)

fig = plt.figure(**figprops) # First sub-plot, the original time series anomaly and inverse wavelet# transform.

ax= plt.subplot(211)
ax.plot(t1, s1, 'k', linewidth=1.5)
ax.plot(t2, s2, 'r', linewidth=1.5)
ax.set_ylabel('RCO', fontdict=font_kk) 
ax.set_xlim(5, 175)
ax.axvspan(SWI_start, SWI_end, facecolor='#005ce6', alpha=0.4) # SW-I
ax.axvspan(SWII_start, SWII_end, facecolor='#006600', alpha=0.4) # SW-II
ax.axvspan(Ploop_start, Ploop_end, facecolor='#ff1a1a', alpha=0.4)  # Ploop

##########################################################################
ex= plt.subplot(212,sharex=ax)
levels = MaxNLocator(nbins=10).tick_values(-1, 1)
cm=ex.contourf(t1,np.log2(cor_period), (WCT), levels=levels)

cbar_ax_1 = fig.add_axes([0.91, 0.11, 0.020, 0.35])
extent3 = [t2.min(), t2.max(), 0, np.log2(max(cor_period))]
ex.contour(t1, np.log2(cor_period), cor_sig, [-99, 1], colors='k', linewidths=2,extent=extent3)
ex.fill(np.concatenate([t2, t2[-1:] + 0, t2[-1:] + 0,t2[:1] - 0, t2[:1] - 0]),
           np.concatenate([np.log2(cor_coi), [1e-9], np.log2(cor_period[-1:]),
           np.log2(cor_period[-1:]), [1e-9]]),'k', alpha=0.3, hatch='x')
ex.quiver(t1[::3], np.log2(cor_period[::3]), u[::3, ::3], v[::3, ::3], units='height',
           angles='uv', pivot='mid', linewidth=0.9, edgecolor='k',
           headwidth=10, headlength=10, headaxislength=5, minshaft=2,
           minlength=5)
ex.set_ylim(1.5, 6.0)
ex.set_xlabel('Res Nos', fontdict=font_kk)
ex.set_ylabel('Period (Arb. units)',fontdict=font_kk)
ex.invert_yaxis()
fig.colorbar(cm, cax=cbar_ax_1)
###############################################################################

fig_file = Protein_name+root_name1+"_"+root_name2+".png"
plt.savefig(fig_file)



fig_file = Protein_name+root_name1+"_"+root_name2+".png"

plt.savefig(fig_file)



figprops = dict(num=2, figsize=(11, 8), dpi=100)

fig = plt.figure(**figprops) # First sub-plot, the original time series anomaly and inverse wavelet# transform.

ax= plt.subplot(111)

yax1 = [0,0]
yax2 = [0, 1.0]
ax.plot(yax1,yax2, color='k')

ax.quiver(0,0, pu,pv, scale=150, facecolor='#005ce6')
ax.quiver(0,0, s1u,s1v, scale=500, facecolor='#006600')
ax.quiver(0,0, s2u,s2v, scale=350, facecolor='#ff1a1a')

ax.set_xlim(-1.0,1.0)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
#ax.set_ylim(-0.01, 0.01)
ax.set_ylim(0.0, 1.0)
ax.set_xlabel('Arb. basis units', fontdict=font_kk)
ax.set_ylabel('Arb. basis units',fontdict=font_kk)
###############################################################################
#ax.text(0.95, 0.15, '$\langle$P_loop-SWI %0.2f'% (Y_ang(pu,pv,s1u,s1v)),
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#3d3c3d', fontsize=18)

#ax.text(0.95, 0.10, '$\langle$P_loop-SWII %0.2f'% (Y_ang(pu,pv,s2u,s2v)),
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#3d3c3d', fontsize=18)

#ax.text(0.95, 0.20, '$\langle$SWI-SWII %0.2f'% (Y_ang(s1u,s1v, s2u,s2v)),
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#3d3c3d', fontsize=18)
ax.text(0.95, 0.15, '$\langle$SWI %0.2f'% (Y_ang(0,1, s1u,s1v)),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='#006600', fontsize=18)
ax.text(0.95, 0.20, '$\langle$SWII %0.2f'% (Y_ang(0,1, s2u,s2v)),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='#ff1a1a', fontsize=18)
ax.text(0.95, 0.25, '$\langle$P-loop %0.2f'% (Y_ang(0,1, pu,pv)),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='#005ce6', fontsize=18)

###############################################################################
#ax.text(0.95, 0.90, 'SWI',
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#006600', fontsize=18)
#
#ax.text(0.95, 0.85, 'SWII',
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#ff1a1a', fontsize=18)
#
#ax.text(0.95, 0.95, 'P-Loop',
#        verticalalignment='bottom', horizontalalignment='right',
#        transform=ax.transAxes,
#        color='#005ce6', fontsize=18)
#
###############################################################################


fig_file = Protein_name+root_name1+"_"+root_name2+"vect"+".png"
plt.savefig(fig_file)




#plt.show()

