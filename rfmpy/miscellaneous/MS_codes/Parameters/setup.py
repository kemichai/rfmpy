import numpy as np


# def Raysum_Params():

#     # SETTING RAYSUM PARAMETERS

#     # Unchanged during the whole execution

#     multiples = 2  # 0 for none, 1 for Moho, 2 for all first-order
#     nSamples = 500  # Traces number of samples
#     syn_delta = 0.05  # Time interval between synthetic samples in seconds
#     syn_fsamp = 1 / syn_delta  # Sampling frequency
#     gauss = 0.6  # Gaussian pulse width (seconds)
#     align = 1  # 0 is none, 1 aligns on primary phase
#     tshift = 10  # Traces shift with respect to 0s
#     rotation = 1  # 0 for N,E,Z - 1 for R,T,Z - 2 for SV,SH,P

#     params = (multiples, nSamples, syn_delta, syn_fsamp, gauss, align, tshift, rotation)

#     return params


def Migration_Params(dx, dy):

    # SETTING GEOMETRY PARAMETERS FOR MIGRATION

    minx = -100 + dx
    maxx = 100 + dx
    pasx = 0.50

    miny = -100 + dy
    maxy = 100 + dy
    pasy = maxy - miny

    minz = -2
    maxz = 100
    pasz = 0.50

    inc = 0.250
    zmax = 100  # !! always <= than maxz

    params = (minx, maxx, pasx, miny, maxy, pasy, minz, maxz, pasz, inc, zmax)

    return params


# def Model_Params(Migration_Params):

# 	# SETTING MODEL PARAMETERS

# 	# They are updated at each misfit evaluation

# 	# Velocity Structure

# 	vS=np.array([3.5,4.5])	# S-wave velocity for upper and lower layer in km/s
# 	vP=np.array([6,7.5]) 		# P-wave velocity for upper and lower layer in km/s
# 	rho=np.array([2700,3100])	# Density for upper and lower layer kg/m3

# 	x1,x2,x3=np.array([46,53,88]) # Model node x-coordinates
# 	z1,z2,z3=np.array([50,5,15])  # Model node z-coordinates

# 	# Seismic interface

# 	minx,maxx,pasx=Migration_Params[:3]
# 	miny,maxy,pasy=Migration_Params[3:6]
# 	minz,maxz,pasz=Migration_Params[6:9]

# 	sx=np.array([minx,x2,x3,maxx])	# Seismic interface nodes | x coord
# 	sz=np.array([z2,z2,z3,z3])		# Seismic interface nodes | z coord

# 	# nbr_segments=len(sx)-1

# 	xModel=np.array([x1,x2,x3])
# 	zModel=np.array([z1,z2,z3])

# 	return xModel,zModel,sx,sz,vS,vP,rho
