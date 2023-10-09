import numpy as np
import os
from matplotlib import pyplot as plt

def between(min, val, max):
    if val < min:
        return min
    elif val > max:
        return max
    else:
        return val

def get_potential(input="LOCPOT"):

    if not os.path.isfile(input):
        print("LOCPOT file not found")
        exit()

    locpot = open(input, 'r')

    #read lattice parameters
    for i in range(2):
        line = locpot.readline()

    lattice = np.zeros((3,3), dtype=float)

    line = locpot.readline()
    lattice[0] = line.split() #Angst
    line = locpot.readline()
    lattice[1] = line.split() #Angst
    line = locpot.readline()
    lattice[2] = line.split() #Angst

    while True: #skip until newline
        line = locpot.readline()
        if not line.strip(): #empty line
            line = locpot.readline()
            break

    #read number of grid points
    grid = np.zeros(3, dtype=int)
    grid[0] = int(line.split()[0])
    grid[1] = int(line.split()[1])
    grid[2] = int(line.split()[2])

    #potential matrix
    potential = np.zeros(grid[0]*grid[1]*grid[2])

    #fill potential matrix
    for i in range(int(grid[0]*grid[1]*grid[2]/5)):
        line = locpot.readline()
        potential[5*i:5*i+5] = np.array(line.split(), dtype=float)

    #last line
    if grid[0]*grid[1]*grid[2]%5 != 0:
        line = locpot.readline()
        potential[-1-len(line.split()):-1] = np.array(line.split(), dtype=float)

    #reshape potential matrix
    potential = np.reshape(potential, (grid[0], grid[1], grid[2]), order='F')
    locpot.close()

    print("Potential data read successfully")
    print("grid points x:" + str(grid[0]) + " y:" + str(grid[1]) + " z:" + str(grid[2]))

    return lattice, grid, potential

def get_charge(input="CHGCAR"):
    if not os.path.isfile(input):
        print("CHGCAR file not found")
        exit()

    locpot = open(input, 'r')

    #read lattice parameters
    for i in range(2):
        line = locpot.readline()

    lattice = np.zeros((3,3), dtype=float)

    line = locpot.readline()
    lattice[0] = line.split() #Angst
    line = locpot.readline()
    lattice[1] = line.split() #Angst
    line = locpot.readline()
    lattice[2] = line.split() #Angst

    while True: #skip until newline
        line = locpot.readline()
        if not line.strip(): #empty line
            line = locpot.readline()
            break

    #read number of grid points
    grid = np.zeros(3, dtype=int)
    grid[0] = int(line.split()[0])
    grid[1] = int(line.split()[1])
    grid[2] = int(line.split()[2])

    #charge matrix
    charge = np.zeros(grid[0]*grid[1]*grid[2])

    #fill charge matrix
    for i in range(int(grid[0]*grid[1]*grid[2]/5)):
        line = locpot.readline()
        charge[5*i:5*i+5] = np.array(line.split(), dtype=float)

    #last line
    if grid[0]*grid[1]*grid[2]%5 != 0:
        line = locpot.readline()
        charge[-1-len(line.split()):-1] = np.array(line.split(), dtype=float)

    #reshape charge matrix
    charge = np.reshape(charge, (grid[0], grid[1], grid[2]), order='F')
    locpot.close()

    print("charge data read successfully")

    return lattice, grid, charge

def draw_stream(axis, slice, name, lattice, grid, potential, charge, broken_stream):
    global V_min, V_max, n_lv

    #lattice constants
    a = np.linalg.norm(lattice[0])
    b = np.linalg.norm(lattice[1])
    c = np.linalg.norm(lattice[2])

    #set up coordinates for plot
    x = np.linspace(0, a, grid[0]+1, endpoint=True)
    y = np.linspace(0, b, grid[1]+1, endpoint=True)
    z = np.linspace(0, c, grid[2]+1, endpoint=True)

    fig1, ax1 = plt.subplots(tight_layout=True)
    
    if axis == 'x':
        x_fig = z
        y_fig = y

        fig1.set_figwidth(between(3, c/4+1, 10))
        fig1.set_figheight(between(2.6, b/4, 10))
        ax1.set_xlim(0,c)
        ax1.set_xlabel("z (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,b)
        ax1.set_ylabel("y (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[1]+1,grid[2]+1))
        potential_fill[:grid[1],:grid[2]] = -potential[slice,:,:]
        potential_fill[grid[1],:grid[2]] = -potential[slice,0,:]
        potential_fill[:grid[1],grid[2]] = -potential[slice,:,0]
        potential_fill[grid[1],grid[2]] = -potential[slice,0,0]
        potential_fill = potential_fill.T

        charge_fill = np.zeros((grid[1]+1,grid[2]+1))
        charge_fill[:grid[1],:grid[2]] = charge[slice,:,:]
        charge_fill[grid[1],:grid[2]] = charge[slice,0,:]
        charge_fill[:grid[1],grid[2]] = charge[slice,:,0]
        charge_fill[grid[1],grid[2]] = charge[slice,0,0]
        charge_fill = charge_fill.T
    elif axis == 'y':
        x_fig = x
        y_fig = z

        fig1.set_figwidth(between(3, a/4+1, 10))
        fig1.set_figheight(between(2.6, c/4, 10))
        ax1.set_xlim(0,a)
        ax1.set_xlabel("x (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,c)
        ax1.set_ylabel("z (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[0]+1,grid[2]+1))
        potential_fill[:grid[0],:grid[2]] = -potential[:,slice,:]
        potential_fill[grid[0],:grid[2]] = -potential[0,slice,:]
        potential_fill[:grid[0],grid[2]] = -potential[:,slice,0]
        potential_fill[grid[0],grid[2]] = -potential[0,slice,0]

        charge_fill = np.zeros((grid[0]+1,grid[2]+1))
        charge_fill[:grid[0],:grid[2]] = charge[:,slice,:]
        charge_fill[grid[0],:grid[2]] = charge[0,slice,:]
        charge_fill[:grid[0],grid[2]] = charge[:,slice,0]
        charge_fill[grid[0],grid[2]] = charge[0,slice,0]
    elif axis == 'z':
        x_fig = x
        y_fig = y

        fig1.set_figwidth(between(3, a/4+1, 10))
        fig1.set_figheight(between(2.6, b/4, 10))
        ax1.set_xlim(0,a)
        ax1.set_xlabel("x (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,b)
        ax1.set_ylabel("y (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[0]+1,grid[1]+1))
        potential_fill[:grid[0],:grid[1]] = -potential[:,:,slice]
        potential_fill[grid[0],:grid[1]] = -potential[0,:,slice]
        potential_fill[:grid[0],grid[1]] = -potential[:,0,slice]
        potential_fill[grid[0],grid[1]] = -potential[0,0,slice]

        charge_fill = np.zeros((grid[0]+1,grid[1]+1))
        charge_fill[:grid[0],:grid[1]] = charge[:,:,slice]
        charge_fill[grid[0],:grid[1]] = charge[0,:,slice]
        charge_fill[:grid[0],grid[1]] = charge[:,0,slice]
        charge_fill[grid[0],grid[1]] = charge[0,0,slice]
    else:
        print("Input axis error")
        return None

    #plot potential contour
    X, Y = np.meshgrid(x_fig, y_fig)
    potential_fill = potential_fill.T
    if n_lv:
        levels = np.linspace(V_min, V_max, n_lv)
        cp1 = ax1.contour(X, Y, potential_fill, levels)
    else:
        cp1 = ax1.contour(X, Y, potential_fill)
        
    cbar1 = fig1.colorbar(cp1)
    cbar1.set_label("potential (V)")
    potential_fill = potential_fill.T

    mask_x = np.array([], dtype=float)
    mask_y = np.array([], dtype=float)
    for i in range(X.shape[0]): #y axis
        for j in range(X.shape[1]): #x axis
            if ((charge_fill > 1)).T[i,j]:
                mask_x = np.concatenate((mask_x, x_fig[j]), axis=None)
                mask_y = np.concatenate((mask_y, y_fig[i]), axis=None)

    ax1.scatter(mask_x,mask_y,s=5,c='black')

    grad = np.gradient(potential_fill, x_fig, y_fig, edge_order=2)
    grad_X = -grad[0]
    grad_Y = -grad[1]

    strength = np.sqrt(grad_X**2 + grad_Y**2)

    lw = 5 * strength / strength.max()
    ax1.streamplot(X, Y, grad_X.T, grad_Y.T, broken_streamlines=broken_stream, density=1, linewidth=1)
    if broken_stream:
        plt.savefig(name + "_stream_broken.png")
    else:
        plt.savefig(name + "_stream.png")

def draw_quiver(axis, slice, name, lattice, grid, potential, charge):
    global V_min, V_max, n_lv

    #lattice constants
    a = np.linalg.norm(lattice[0])
    b = np.linalg.norm(lattice[1])
    c = np.linalg.norm(lattice[2])

    #set up coordinates for plot
    x = np.linspace(0, a, grid[0]+1, endpoint=True)
    y = np.linspace(0, b, grid[1]+1, endpoint=True)
    z = np.linspace(0, c, grid[2]+1, endpoint=True)

    fig1, ax1 = plt.subplots(tight_layout=True)

    if axis == 'x':
        x_fig = z
        y_fig = y

        fig1.set_figwidth(between(3, c/4+1, 10))
        fig1.set_figheight(between(2.6, b/4, 10))
        ax1.set_xlim(0,c)
        ax1.set_xlabel("z (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,b)
        ax1.set_ylabel("y (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[1]+1,grid[2]+1))
        potential_fill[:grid[1],:grid[2]] = -potential[slice,:,:]
        potential_fill[grid[1],:grid[2]] = -potential[slice,0,:]
        potential_fill[:grid[1],grid[2]] = -potential[slice,:,0]
        potential_fill[grid[1],grid[2]] = -potential[slice,0,0]
        potential_fill = potential_fill.T

        charge_fill = np.zeros((grid[1]+1,grid[2]+1))
        charge_fill[:grid[1],:grid[2]] = charge[slice,:,:]
        charge_fill[grid[1],:grid[2]] = charge[slice,0,:]
        charge_fill[:grid[1],grid[2]] = charge[slice,:,0]
        charge_fill[grid[1],grid[2]] = charge[slice,0,0]
        charge_fill = charge_fill.T
    elif axis == 'y':
        x_fig = x
        y_fig = z

        fig1.set_figwidth(between(3, a/4+1, 10))
        fig1.set_figheight(between(2.6, c/4, 10))
        ax1.set_xlim(0,a)
        ax1.set_xlabel("x (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,c)
        ax1.set_ylabel("z (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[0]+1,grid[2]+1))
        potential_fill[:grid[0],:grid[2]] = -potential[:,slice,:]
        potential_fill[grid[0],:grid[2]] = -potential[0,slice,:]
        potential_fill[:grid[0],grid[2]] = -potential[:,slice,0]
        potential_fill[grid[0],grid[2]] = -potential[0,slice,0]

        charge_fill = np.zeros((grid[0]+1,grid[2]+1))
        charge_fill[:grid[0],:grid[2]] = charge[:,slice,:]
        charge_fill[grid[0],:grid[2]] = charge[0,slice,:]
        charge_fill[:grid[0],grid[2]] = charge[:,slice,0]
        charge_fill[grid[0],grid[2]] = charge[0,slice,0]
    elif axis == 'z':
        x_fig = x
        y_fig = y

        fig1.set_figwidth(between(3, a/4+1, 10))
        fig1.set_figheight(between(2.6, b/4, 10))
        ax1.set_xlim(0,a)
        ax1.set_xlabel("x (" + 	u"\u212B" + ")")
        ax1.set_ylim(0,b)
        ax1.set_ylabel("y (" + 	u"\u212B" + ")")

        #reshape based on repeating lattice
        potential_fill = np.zeros((grid[0]+1,grid[1]+1))
        potential_fill[:grid[0],:grid[1]] = -potential[:,:,slice]
        potential_fill[grid[0],:grid[1]] = -potential[0,:,slice]
        potential_fill[:grid[0],grid[1]] = -potential[:,0,slice]
        potential_fill[grid[0],grid[1]] = -potential[0,0,slice]

        charge_fill = np.zeros((grid[0]+1,grid[1]+1))
        charge_fill[:grid[0],:grid[1]] = charge[:,:,slice]
        charge_fill[grid[0],:grid[1]] = charge[0,:,slice]
        charge_fill[:grid[0],grid[1]] = charge[:,0,slice]
        charge_fill[grid[0],grid[1]] = charge[0,0,slice]
    else:
        print("Input axis error")
        return None

    print("\nmax potential:" + "{:g}".format(np.max(potential_fill)))
    print("min potential:" + "{:g}".format(np.min(potential_fill)))
    print()

    V_min = float(input("Minimum potential (optional): ") or "NaN")

    if not np.isnan(V_min):
        V_max = float(input("Maximum potential: ") or "NaN")
    else: 
        V_max = float("NaN")

    if not np.isnan(V_max):
        n_lv = int(input("Number of potential levels: ") or 0)
    else: 
        n_lv = 0

    #plot potential contour
    X, Y = np.meshgrid(x_fig, y_fig)
    potential_fill = potential_fill.T
    if n_lv:
        levels = np.linspace(V_min, V_max, n_lv)
        cp1 = ax1.contour(X, Y, potential_fill, levels)
    else:
        cp1 = ax1.contour(X, Y, potential_fill)
        
    cbar1 = fig1.colorbar(cp1)
    cbar1.set_label("potential (V)")
    potential_fill = potential_fill.T

    grad = np.gradient(potential_fill, x_fig, y_fig, edge_order=2)
    grad_X = -grad[0]
    grad_Y = -grad[1]

    #mask areas with high charges
    grad_X[(charge_fill > 1)] = 0
    grad_Y[(charge_fill > 1)] = 0

    #mask areas at the dipole correction
    if axis == 'x': 
        grad_X[(X.T < 0 + 0.3) | (X.T > c - 0.3)] = 0
        grad_Y[(X.T < 0 + 0.3) | (X.T > c - 0.3)] = 0
    elif axis == 'y': 
        grad_X[(Y.T < 0 + 0.3) | (Y.T > c - 0.3)] = 0
        grad_Y[(Y.T < 0 + 0.3) | (Y.T > c - 0.3)] = 0

    mask_x = np.array([], dtype=float)
    mask_y = np.array([], dtype=float)
    for i in range(X.shape[0]): #y axis
        for j in range(X.shape[1]): #x axis
            if ((charge_fill > 1)).T[i,j]:
                mask_x = np.concatenate((mask_x, x_fig[j]), axis=None)
                mask_y = np.concatenate((mask_y, y_fig[i]), axis=None)

    ax1.scatter(mask_x,mask_y,s=6,c='black')

    step=8
    ax1.quiver(X[::step,::step], Y[::step,::step], grad_X.T[::step,::step], grad_Y.T[::step,::step], width=0.002, scale=8)    
    plt.savefig(name + "_contour.png")


    # fig, ax = plt.subplots(tight_layout=True)
    # cp = ax.imshow(potential_fill, vmin=V_min, vmax=V_max)
    # cbar = fig.colorbar(cp)
    # cbar.set_label("potential (V)")
    # X, Y = np.meshgrid(np.arange(grid[0]+1), np.arange(grid[1]+1))
    # plt.quiver(X, Y, grad_X.T, -grad_Y.T) #imshow plots from top to bottom
    # plt.savefig(name + "_heatmap.png")

#start here!

lattice, grid, potential = get_potential()
_, _, charge = get_charge()

axis = input("\ninput axis to slice (x, y or z): ")
slice = int(input("input plane index to slice (0 to n-1): "))
name = input("name of the plane: ")

draw_quiver(axis, slice, name, lattice, grid, potential, charge)
draw_stream(axis, slice, name, lattice, grid, potential, charge, broken_stream=True)
#draw_stream(axis, slice, name, lattice, grid, potential, charge, False)
plt.show()