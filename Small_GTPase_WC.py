#!/usr/bin/env python
import matplotlib.pyplot as plt
import csv
import os
import numpy as np
import scipy as sp
import sys
import pycwt as wavelet
from pycwt.helpers import find
from matplotlib.image import NonUniformImage
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D


###################################################################################################################################
#PLEASE NOTE: YOU HAVE TO EDIT THE PATH OF WCT SIGNIFICANCE FILE PATH IN " site-packages/pycwt/wavelets.py" (Line numbers 568 & 632).
#IT SHOULD BE SAME AS IN THE CURRENT FILE, LINE NUMEBRS 132, 135 AND 215
################################################################################################################################## 

while True:
	sys.stdout.write('Have you edited this script for start and end residue numbers of P-loop, SWI, SWII regions?\n')
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
SWI_end = 37

SWII_start = 59
SWII_end = 76

Ploop_start = 12
Ploop_end = 18
###################################################################################################

SMALL_SIZE = 14.5
MEDIUM_SIZE = 16.5
BIGGER_SIZE = 18.5

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



#scg.set_default_wavelet('cmor1-1.5')
#scg.set_default_wavelet('morl')


a = []
a1 = []
a3 = []
a4 = []


root_name1 = sys.argv[1]

file_name1=root_name1+".csv"


if os.path.exists(file_name1):
	dum=1
else:
  print("The file %s does not exist" % file_name1)
  exit()

with open(file_name1,'rU') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
                a.append(int(row[0]))
                a1.append(float(row[1]))


t1 = np.asarray(a) 
s1 = np.asarray(a1) 

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

with open(file_name2,'rU') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
                a3.append(int(row[0]))
                a4.append(float(row[1]))


t2 = np.asarray(a3) 
s2 = np.asarray(a4) 


if (np.array_equal(t1,t2)):
	print("The two signals are of identical lengths")
else:
	comb_name = file_name1+","+file_name2
	print("CSV files have diffrent lenths %s " % comb_name)
	exit()

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
s0 = -1  # 2 * dt                   # Starting scale, here 6 months
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

#################################################################################################

os.rename('/Users/admin/.cache/pycwt/kl.txt', sim_file )

U1 = u[SWI_start:SWI_start+1:1, SWI_start:SWI_start+1:1]
for i in range(SWI_start,SWI_end):
	U1 = U1 + u[i:i+1:1, i:i+1:1]
	print(U1)

V1 = v[SWI_start:SWI_start+1:1, SWI_start:SWI_start+1:1]
for i in range(SWI_start,SWI_end):
	V1 = V1 + v[i:i+1:1, i:i+1:1]

U2 = u[SWII_start:SWII_start+1:1, SWII_start:SWII_start+1:1]
for i in range(SWII_start,SWII_end):
	U2 = U2 + u[i:i+1:1, i:i+1:1]

V2 = v[SWII_start:SWII_start+1:1, SWII_start:SWII_start+1:1]
for i in range(SWII_start,SWII_end):
	V2 = V2 + v[i:i+1:1, i:i+1:1]

U3 = u[Ploop_start:Ploop_start+1:1, Ploop_start:Ploop_start+1:1]
for i in range(Ploop_start, Ploop_end):
	U3 = U3 + u[i:i+1:1, i:i+1:1]

V3 = v[Ploop_start:Ploop_start+1:1, Ploop_start:Ploop_start+1:1]
for i in range(Ploop_start, Ploop_end):
	V3 = V3 + v[i:i+1:1, i:i+1:1]

yax1 = [0,0]
yax2 = [0,0.01]
yax = np.array([yax1, yax2])
l1_a = np.array([0, U1])
l1_b = np.array([0, V1])
l2_a = np.array([0, U2])
l2_b = np.array([0, V2])

line_1 = Line2D(l1_a, l1_b, linewidth=1, linestyle = "-", color="green")
line_2 = Line2D(l2_a, l2_b, linewidth=1, linestyle = "-", color="green")

#print((cor_power))



########################Calculate the angle between vectors########################


import math

def GetAngle(x1, y1, x2,y2):
    inner_product = x1*x2 + y1*y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    return math.degrees(math.acos(inner_product/(len1*len2)))

print(" Angle between SWI and SW2 %.2f" % GetAngle(U1, V1, U2, V2))
print(" Angle between SWI and Ploop %.2f" % GetAngle(U1, V1, U3, V3))
print(" Angle between SWII and Ploop %.2f" % GetAngle(U2, V2, U3, V3))


######################################Plot the Signal(RCO)##############################
plt.ioff()
figprops = dict(num=1, figsize=(11, 8), dpi=300)

fig = plt.figure(**figprops) # First sub-plot, the original time series anomaly and inverse wavelet# transform.

ax= plt.subplot(211)
ax.plot(t1, s1, 'k', linewidth=1.5)
ax.plot(t2, s2, 'r', linewidth=1.5)
ax.set_ylabel('RCO') 
ax.set_xlabel('Res Nos')
ax.axvspan(SWI_start, SWI_end, facecolor='#005ce6', alpha=0.4) # SW-I
ax.axvspan(SWII_start, SWII_end, facecolor='#006600', alpha=0.4) # SW-II
ax.axvspan(Ploop_start, Ploop_end, facecolor='#ff1a1a', alpha=0.4)  # Ploop

extent_corr = [t1.min(), t1.max(), 0, np.log2(max(cor_period))]
ex= plt.subplot(212,sharex=ax)

##########################################################################
levels = MaxNLocator(nbins=10).tick_values(-1, 1)
cm=ex.contourf(t1,
                  #np.log2(cor_period), np.log2(cor_power), levels=levels,
                  np.log2(cor_period), (WCT), levels=levels,
                  #(cor_period), (cor_power), levels=levels,
                  )

###############################################################################
cbar_ax_1 = fig.add_axes([0.91, 0.11, 0.020, 0.35])

extent3 = [t2.min(), t2.max(), 0, np.log2(max(cor_period))]
ex.contour(t2, np.log2(cor_period), cor_sig, [-99, 1], colors='k', linewidths=2,extent=extent3)
ex.fill(np.concatenate([t2, t2[-1:] + 0, t2[-1:] + 0,t2[:1] - 0, t2[:1] - 0]),np.concatenate([np.log2(cor_coi), [1e-9], np.log2(cor_period[-1:]),np.log2(cor_period[-1:]), [1e-9]]),'k', alpha=0.3, hatch='x')
ex.quiver(t1[::3], np.log2(cor_period[::3]), u[::3, ::3], v[::3, ::3], units='height',
           angles='uv', pivot='mid', linewidth=0.9, edgecolor='k',
           headwidth=10, headlength=10, headaxislength=5, minshaft=2,
           minlength=5)
ex.set_ylim(1.5, 6.0)
ex.set_xlabel('Res Nos')
ex.invert_yaxis()
fig.colorbar(cm, cax=cbar_ax_1)

#print(WCT)


fig_file = Protein_name+root_name1+"_"+root_name2+".png"

plt.savefig(fig_file)



figprops = dict(num=2, figsize=(11, 8), dpi=400)

fig = plt.figure(**figprops) # First sub-plot, the original time series anomaly and inverse wavelet# transform.

ax= plt.subplot(111)

ax.plot(yax1,yax2, color='k')
ax.quiver(0,0, U1, V1, scale=20, facecolor='#005ce6') 
ax.quiver(0,0, U2, V2, scale=20, facecolor='#006600') 
ax.quiver(0,0, U3, V3, scale=20, facecolor='#ff1a1a') 
ax.set_xlim(-0.01,0.01)
#ax.set_ylim(-0.01, 0.01)
ax.set_ylim(0.0, 0.01)
ax.text(0.95, 0.01, '$\langle$SWI-SWII %0.2f'% GetAngle(U1, V1, U2, V2),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=15)

ax.text(0.95, 0.04, '$\langle$SWI-Ploop %0.2f'% GetAngle(U1, V1, U3, V3),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=15)

ax.text(0.95, 0.07, '$\langle$SWII-Ploop %0.2f'% GetAngle(U2, V2, U3, V3),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=15)

fig_file = Protein_name+"_"+root_name1+"_"+root_name2+"vect"+".png"
fig_file = Protein_name+"_"+root_name1+"_"+root_name2+"vect"+".png"
fig_file = Protein_name+"_"+root_name1+"_"+root_name2+"vect"+".png"

plt.tight_layout()
plt.savefig(fig_file)


#plt.show()

