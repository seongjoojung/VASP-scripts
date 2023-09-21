"""
NOTE: only implemented for orthogonal systems  

Requires LCALCEPS = T tag in VASP INCAR
Requires CONTCAR and OUTCAR in directory
"""   
import numpy as np
import sys
import os
import re

def get_energy(outcar="OUTCAR"):
    if not os.path.isfile(outcar):
        print("OUTCAR file not found")
        return None
    
    txt = open(outcar).read()
    energy = re.findall(r'energy\(sigma->0\) =\s*[-+]?(?:\d*\.*\d+)', txt)[-1].split()[-1]

    return float(energy)

def get_dipoles(outcar="OUTCAR"):
    if not os.path.isfile(outcar):
        print("OUTCAR file not found")
        return None
    
    txt = open(outcar).read()
    dipole_ion = re.findall(r'p\[ion\]=\(\s*[-+]?(?:\d*\.*\d+)\s*[-+]?(?:\d*\.*\d+)\s*[-+]?(?:\d*\.*\d+)\s*\)', txt)[-1].split()[1:4]
    dipole_ion = np.array(dipole_ion, dtype=float)
    dipole_elec = re.findall(r'p\[elc\]=\(\s*[-+]?(?:\d*\.*\d+)\s*[-+]?(?:\d*\.*\d+)\s*[-+]?(?:\d*\.*\d+)\s*\)', txt)[-1].split()[1:4]
    dipole_elec = np.array(dipole_elec, dtype=float)

    return dipole_ion, dipole_elec

def get_lattice(contcar="CONTCAR"):
    if not os.path.isfile(contcar):
        print("CONTCAR file not found")
        return None

    #read CONTCAR
    contcar = open('CONTCAR', 'r')

    #skip first 2 lines
    for i in range(3):
        line = contcar.readline()

    #lattice constants
    a = float(line.split()[0]) #Angst
    line = contcar.readline()
    b = float(line.split()[1]) #Angst
    line = contcar.readline()
    c = float(line.split()[2]) #Angst
    contcar.close()

    lattice = np.array([a,b,c], dtype=float)

    return lattice

#start here
n_arg = len(sys.argv)

if n_arg < 4: 
    print("Provide P quantum offset for Berry phase polarization calculation")
elif n_arg == 4: #non-polar ionic dipole moment provided
    print("P quantum offset read: " + sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3])
    P_quantum_offset = np.array([sys.argv[1], sys.argv[2], sys.argv[3]], dtype=float)
else: 
    print("Provide P quantum offset for Berry phase polarization calculation")

#read files
lattice = get_lattice()
a = lattice[0]
b = lattice[1]
c = lattice[2]
dipole_ion, dipole_elec = get_dipoles()

#elementary charge
e = 1.602*10**-19 #C

#calculate Berry phase polarization
np.set_printoptions(precision=5,suppress=True)
print("ionic dipole:")
print(dipole_ion)
print("electronic dipole:")
print(dipole_elec)

P_quantum = lattice*(e*10**20)/(a*b*c)
print("Polarization quantum (C/m^2):")
print(P_quantum)

#dipole moment is in units of electrons*Angst, opposite sign!
Berry_P = -(dipole_ion + dipole_elec)*(e*10**20)/(a*b*c)

for i in range(3):
    while(Berry_P[i] < -1e-4):
        Berry_P[i] += P_quantum[i]

    while(Berry_P[i] > P_quantum[i] - 1e-4):
        Berry_P[i] -= P_quantum[i]

#add offset
Berry_P += P_quantum_offset*P_quantum

print("\nBerry phase polarization (C/m^2):")
print(Berry_P)

print("Other possible values: ")
print(Berry_P - P_quantum)
print(Berry_P + P_quantum)

np.savetxt("Berry_P.dat", Berry_P, fmt='%.5f')
exit
