# Code to take original and "spherified" high-res dipole positions and visualise them in 3D.
#
# INSTRUCTIONS FOR USE:
#
# - See 'spherify' for instructions, as this should only be run as part of this package.
#
# The code works on the mac it was designed on, but has not been tested on any other computers - please e-mail
# the following e-mail address if you have any questions or issues: matthew.g.lodge@gmail.com
#
# TROUBLESHOOTING
#
# This exact version is designed to produce two windows (side-by-side) on mac computers; however, this will not
# work if running on windows. Simply comment out any lines that produce errors as a workaround, e.g:
#
# - matplotlib.use("TkAgg") # use a backend to allow the current_fig_manager section to position the windows
# - plt.get_current_fig_manager().window.wm_geometry("+750+100")
#
# and then uncomment the lines under "windows version".
# 
# Matt Lodge 04/08/22

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import pandas as pd
matplotlib.use("TkAgg") # use a backend to allow the current_fig_manager section to position the windows
import matplotlib.pyplot as plt

print("\n\n ------------------------------------------------------------------------------------------------")
print("     Welcome to S.T.A.G (Simulated Three-dimensional Aerosol Geometries!): Spherify Edition")
print(" ------------------------------------------------------------------------------------------------")

# ------- ORIGINAL POSITIONS -------

print("\n ORIGINAL POSITIONS:\n")

# Import dipole positions
original = pd.read_csv('original.txt', header=None, names=['X', 'Y', 'Z'])

print(" ",len(original)-1," dipoles imported successfully from original image.")


# The final row gives information about the lattice dimensions and the number of dipoles - store this info before moving on
STAG_lattice_dim=original['X'][len(original)-1] # store the lattice dimension value
N=original['Y'][len(original)-1] # store the number of dipoles N from our program


# create a 3D grid composed entirely of zeroes, with our lattice dimensions
grid= np.zeros((STAG_lattice_dim,STAG_lattice_dim,STAG_lattice_dim),dtype=int) # create a grid composed of zeroes
print(" ",STAG_lattice_dim,"x",STAG_lattice_dim,"x",STAG_lattice_dim," grid created.\n")


# Now that we know the number of dipoles "N", scan through all rows of the imported integer matrix "dipoles" and change any points in our zero-array "grid" from 0->1 wherever there are dipole positions (at any dipole coordinate)
for i in range (N):
    grid[original['X'][i]][original['Y'][i]][original['Z'][i]]=1  


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#mac version
plt.get_current_fig_manager().window.wm_geometry("+50+100")

#windows verison:
#mngr = plt.get_current_fig_manager()
#mngr.window.setGeometry(50,100,900,900)

# set axis titles
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# plot voxels
ax.voxels(grid, edgecolor="k")

print(" Original image loaded.\n\n")



ax.set_axis_off() # optional: removes axes and grey area around shape

# uncomment this line to print images seperately:   plt.show()

# ------- HIGH-RESOLUTION POSITIONS -------

print(" SPHERIFIED POSITIONS:\n")

# Import dipole positions
high_res = pd.read_csv('high_res.txt', header=None, names=['X', 'Y', 'Z'])

print(" ",len(high_res)-1," dipoles imported successfully from spherified image.")


# The final row gives information about the lattice dimensions and the number of dipoles - store this info before moving on
STAG_lattice_dim=high_res['X'][len(high_res)-1] # store the lattice dimension value
N=high_res['Y'][len(high_res)-1] # store the number of dipoles N from our program

# create a 3D grid composed entirely of zeroes, with our lattice dimensions
new_grid= np.zeros((STAG_lattice_dim,STAG_lattice_dim,STAG_lattice_dim),dtype=int) # create a grid composed of zeroes
print(" ",STAG_lattice_dim,"x",STAG_lattice_dim,"x",STAG_lattice_dim," grid created.\n")


# Now that we know the number of dipoles "N", scan through all rows of the imported integer matrix "dipoles" and change any points in our zero-array "grid" from 0->1 wherever there are dipole positions (at any dipole coordinate)
for i in range (N):
    new_grid[high_res['X'][i]][high_res['Y'][i]][high_res['Z'][i]]=1  


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#mac version
plt.get_current_fig_manager().window.wm_geometry("+750+100")

#windows verison:
#mngr = plt.get_current_fig_manager()
#mngr.window.setGeometry(950,100,900,900)

# set axis titles
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# plot voxels
ax.voxels(new_grid, edgecolor="k")

print(" Spherified image loaded.\n\n")

plt.show()

