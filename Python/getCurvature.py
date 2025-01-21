import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

df = pd.read_csv('Matlab/code/SingleSectionCR_L0.60_N100_r0.0100_Tt10.50_Tt20.00.csv',header=None)

df2 = pd.read_csv('Matlab/code/SingleSectionCR_L0.60_N100_r0.0100_Tt10.40_Tt20.00.csv',header=None)
df3 = pd.read_csv('Matlab/code/SingleSectionCR_L0.60_N100_r0.0100_Tt10.30_Tt20.00.csv',header=None)
dfArr = [df, df2, df3]

def getCurvature(df):
    """
    Get the curvature of the track.
    """
    x = []
    z = []

    for i, row in df.iterrows():
        x.append(row[0])
        z.append(row[2])

    dx_dt = np.gradient(np.array(x))
    dz_dt = np.gradient(np.array(z))
    velocity = np.array([[dx_dt[i], dz_dt[i]] for i in range(dx_dt.size)])

    ds_dt = np.sqrt(dx_dt * dx_dt + dz_dt * dz_dt)
    
    tangent = np.array([1/ds_dt] * 2).transpose() * velocity
    
    tangent_x = tangent[:,0]
    tangent_z = tangent[:,1]
    
    deriv_tangent_x = np.gradient(tangent_x)
    deriv_tangent_z = np.gradient(tangent_z)
    
    dT_dt = np.array([ [deriv_tangent_x[i], deriv_tangent_z[i]] for i in range
    (deriv_tangent_x.size)])
    
    length_dT_dt = np.sqrt(deriv_tangent_x * deriv_tangent_x + deriv_tangent_z * deriv_tangent_z)
    
    normal = np.array([1/length_dT_dt] * 2).transpose() * dT_dt

    d2s_dt2 = np.gradient(ds_dt)
    d2x_dt2 = np.gradient(dx_dt)
    d2z_dt2 = np.gradient(dz_dt)
    curvature = (d2x_dt2 * dz_dt - dx_dt * d2z_dt2) / ((dx_dt*dx_dt + dz_dt*dz_dt)**1.5)
    return curvature


curvatureArr = []
for df in dfArr:
    curvatureArr.append(getCurvature(df))
print(curvatureArr)
for array in curvatureArr:
    print(array)
    plt.plot(np.linspace(0,0.6,100-5), -array[5:])
plt.xlabel('Position along the track')
plt.ylabel('Curvature')
plt.title('Curvature along the track')
plt.grid(True)
plt.show()
plt.plot(x,z)
plt.xlabel('x')
plt.ylabel('z')
plt.axis([0,1,-1,1])
plt.grid(True)
plt.show()