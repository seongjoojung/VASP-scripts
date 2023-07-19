import numpy as np
import sys
import warnings

# Ignore complex to float warning
warnings.filterwarnings('ignore')

# small values to ignore
n = len(sys.argv)
if n == 4:
    zero_NMC = float(sys.argv[1])
    zero_MOS = float(sys.argv[2])
    zero_ICP = float(sys.argv[3])
else:
    zero_NMC = float(input("zero criterion for Normal Mode Charges: "))
    zero_MOS = float(input("zero criterion for Mode Oscillator Strength: "))
    zero_ICP = float(input("zero criterion for Ionic Contribution to Permittivity: "))

#############################################
#     Read data from CONTCAR and OUTCAR     #
#############################################

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
volume = a*b*c*10**-30     #m^-3

#read elements from next line
line = contcar.readline()
elements = line.split()

#read number of atoms for each element from next line
line = contcar.readline()
element_count = np.array(line.split(), dtype=int)
n_ion = np.sum(element_count)

contcar.close()

#Read OUTCAR
outcar = open('OUTCAR', 'r')

#skip until you find 'POMASS'
while True:
    line = outcar.readline()
    if len(line.split()) == 0:
        continue
    if line.split()[0] == 'POMASS':
        break

atomic_weight = np.array([]) #unitless
for i in range(len(elements)):
    for j in range(element_count[i]):
        for k in range(3):
            atomic_weight = np.append(atomic_weight, float(line.split()[2+i]))

#infinite frequency dielectric tensor
eps_inf = np.zeros((3,3))

#skip until you find 'MACROSCOPIC STATIC DIELECTRIC TENSOR' twice
while True:
    line = outcar.readline()
    if len(line.split()) == 0:
        continue
    if line.split()[0] == 'MACROSCOPIC':
        outcar.readline()
        break

while True:
    line = outcar.readline()
    if len(line.split()) == 0:
        continue
    if line.split()[0] == 'MACROSCOPIC':
        outcar.readline()
        break

for i in range(3):
    line = outcar.readline()

    for j in range(3):
        eps_inf[i,j] = float(line.split()[j])

np.savetxt("permittivity_inf.dat", eps_inf, fmt='%.5g', delimiter='\t')

#Born effective charge tensor
BEC = np.zeros((n_ion, 3, 3))

#skip until you find 'BORN EFFECTIVE CHARGES'
while True:
    line = outcar.readline()
    if len(line.split()) == 0:
        continue
    if line.split()[0] == 'BORN':
        outcar.readline()
        break

for n in range(n_ion):
    outcar.readline()
    for i in range(3):
        line = outcar.readline()
        for j in range(3):
            BEC[n,i,j] = float(line.split()[1+j])

with open('BEC.dat', 'w') as f:
    for n in range(n_ion):
        f.write('ion: ' + str(n) + '\n')
        np.savetxt(f, BEC[n], fmt='%.5g', delimiter="\t")

np.savetxt("BEC_z.dat", BEC[:,2,2], fmt='%.5g', delimiter='\t')

eigenvalues = np.zeros(3*n_ion, dtype=complex) #2PiTHz
eigenvectors = np.zeros((3*n_ion,3*n_ion))     #each row is eigenvector, not column

#skip until you find 'Eigenvectors' again
while True:
    line = outcar.readline()
    if len(line.split()) == 0:
        continue
    if line.split()[0] == 'Eigenvectors':
        outcar.readline()
        outcar.readline()
        outcar.readline()
        break

for i in range(3*n_ion):
    line = outcar.readline()

    if line.split()[1] == 'f':
        eigenvalues[i] = float(line.split()[5])
    else:
        eigenvalues[i] = float(line.split()[4])*1j

    outcar.readline()

    for j in range(n_ion):
        line = outcar.readline()
        eigenvectors[i,3*j+0] = float(line.split()[3])
        eigenvectors[i,3*j+1] = float(line.split()[4])
        eigenvectors[i,3*j+2] = float(line.split()[5])

    outcar.readline()

print("\nEigenvalues (2PiTHz):")
print(eigenvalues)

np.savetxt("eigenvalues.dat", eigenvalues, fmt='%.5g', delimiter='\t')
np.savetxt("eigenvectors.dat", eigenvectors, fmt='%.5g', delimiter='\t')

eigenvalues = eigenvalues*10**12 #2PiHz

#############################################
#               Calculate stuff             #
#############################################

e_charge = 1.602177*10**-19  #C
vac_per = 8.854*10**-12      #F/m
proton_mass = 1.673*10**-27  #kg

U = eigenvectors/(np.sqrt(proton_mass*atomic_weight)) #kg^-1/2

N = np.zeros((3*n_ion, 3*n_ion)) #normal mode, unitless
for i in range(3*n_ion):
    N[i] = U[i]/np.linalg.norm(U[i])

NMC = np.zeros((3*n_ion, 3)) #normal mode charge

#reshape BEC tensor to 3x15 matrix
BEC_reshape = np.zeros((3, 3*n_ion))
for i in range(n_ion):
    BEC_reshape[:,3*i:3*i+3] = BEC[i]

for i in range(3*n_ion):
    NMC[i] = BEC_reshape @ N[i]

#zero out the error
NMC[np.abs(NMC) < zero_NMC] = 0
np.set_printoptions(precision=3)
print("\nNormal Mode Charges:")
print(NMC)

MOS = np.zeros((3*n_ion, 3, 3)) #Mode oscilator strength, C^2/kg

for m in range(3*n_ion):
    a = U[m] @ (BEC_reshape.T*e_charge)
    a = np.reshape(a, (1,3))
    MOS[m] = a.T @ a

#zero out the error
MOS[np.abs(MOS) < zero_MOS] = 0

print("\nMode Oscillator Strength (C^2/kg):")
print(MOS)

ICP = np.zeros((3,3)) #ionic contribution to permittivity

for m in range(3*n_ion):
    ICP += 1/(volume*vac_per)*MOS[m]/(float(eigenvalues[m]**2))

#zero out the error
ICP[np.abs(ICP) < zero_ICP] = 0

print("\nIonic Contribution to Dielectric Constant:", end='\n')
print(ICP, end='\n')

np.savetxt("permittivity_ionic.dat", ICP, fmt='%.5g', delimiter='\t')