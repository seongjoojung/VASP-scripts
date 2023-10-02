"""
noncollinear PDOS 
using ElementTree to import xml data
"""

import xml.etree.ElementTree as ET
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#read inputs
ion_text = input("List of ions (whitespace-separated, index starts from 0):")

try:
    ions = np.array(ion_text.split(), dtype=int)
except ValueError:
    raise ValueError('Input must be integers')

try:
    x_low = float(input("energy lower limit (default = -10): ") or -10) 
except ValueError:
    raise ValueError('Input must be float')

try:
    x_upp = float(input("energy upper limit (default = 5): ") or 5) 
except ValueError:
    raise ValueError('Input must be float')
    
try:    
    y_limit = float(input("DOS plot limit (default = 1.2*max): ") or "NaN") 
except ValueError:
    raise ValueError('Input must be float')

#number of ions
N = ions.shape[0]

tree = ET.parse('vasprun.xml')
root = tree.getroot()
dos = list(root.iter("dos"))[0] #root.iter returns iteratable, list()[0] returns element

efermi = float(dos[0].text.split()[0])

nrow = len(dos[2][0][13][0][0])
ncol = len(dos[2][0][13][0][0][0].text.split())
DOS = np.zeros((N, nrow, ncol))

TDOS = np.zeros((nrow, 3))

#0:x 1:s, 2:py, 3:pz, 4:px, 5:dxy, 6:dyz, 7:dz2, 8:dxz, 9:dx2-y2
#atom number dos[2][0][13][j][0][i], j starts from 0
for i in range(nrow):
    TDOS[i] = np.array(dos[1][0][5][0][i].text.split(), dtype=float)
    for n, j in enumerate(ions): #n: 0-N atom index, j: atom number
        DOS[n][i] = np.array(dos[2][0][13][j][0][i].text.split(), dtype=float) 

#split x out
x = DOS[0][:,0]
x -= efermi #set fermi level to 0
DOS = DOS[:,:,1:]
#0:s, 1:py, 2:pz, 3:px, 4:dxy, 5:dyz, 6:dz2, 7:dxz, 8:dx2-y2

for n, j in enumerate(ions): #n: 0-N atom index, j: atom number

    cm = 1/2.54
    fig = plt.figure(figsize=(18*cm,5*cm), tight_layout=False)
    s_DOS = DOS[n][:, 0]
    p_DOS = DOS[n][:, 1]+DOS[n][:, 2]+DOS[n][:, 3]
    d_DOS = DOS[n][:, 4]+DOS[n][:, 5]+DOS[n][:, 6]+DOS[n][:, 7]+DOS[n][:, 8]
    ax = plt.axes()
    plt.plot(x, s_DOS, label="$s$")
    plt.plot(x, p_DOS, label="$p$")
    plt.plot(x, d_DOS, label="$d$")
    plt.xlim([x_low,x_upp])

    #find maximum value of DOS
    DOS_max = max(np.max(s_DOS[(x>x_low) & (x<x_upp)]), np.max(p_DOS[(x>x_low) & (x<x_upp)]), np.max(d_DOS[(x>x_low) & (x<x_upp)]))
    
    if np.isnan(y_limit):
        plt.ylim([0,DOS_max*1.2])
    else:
        plt.ylim([0,y_limit])

    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(axis="x", which='major', direction='in', length=6, labelsize=11)
    ax.tick_params(axis="x", which='minor', direction='in', length=4)
    
    plt.vlines(0, ax.get_ylim()[0],ax.get_ylim()[1], 'gray', linestyles='dotted')
 
    plt.xlabel("$E-E_f$ (eV)", fontsize=13) 
    
    plt.ylabel("DOS (a.u.)", fontsize=13)   
    plt.legend(prop={'size': 11}, frameon=False)
    plt.savefig("Ion_" + str(j) + ".png" , transparent=True)
    plt.show()

