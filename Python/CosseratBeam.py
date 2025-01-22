import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class CosseratBeam():
    def __init__(self, filelocation):
        self.df = pd.read_csv(filelocation,header=None)
        self.x = []
        self.z = []
        for i, row in self.df.iterrows():
            self.x.append(row[0])
            self.z.append(row[2])

    def getCurvature(self):
        """
        Get the curvature of the track.
        """
        x = []
        z = []

        for i, row in self.df.iterrows():
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
        
        dT_dt = np.array([ [deriv_tangent_x[i], deriv_tangent_z[i]] for i in range(deriv_tangent_x.size)])
        
        length_dT_dt = np.sqrt(deriv_tangent_x * deriv_tangent_x + deriv_tangent_z * deriv_tangent_z)
        
        normal = np.array([1/length_dT_dt] * 2).transpose() * dT_dt

        d2s_dt2 = np.gradient(ds_dt)
        d2x_dt2 = np.gradient(dx_dt)
        d2z_dt2 = np.gradient(dz_dt)
        curvature = (d2x_dt2 * dz_dt - dx_dt * d2z_dt2) / ((dx_dt*dx_dt + dz_dt*dz_dt)**1.5)
        return curvature
    
    def getCurvatureNew(self):
        """
        Get the curvature of the track.
        """
        curvature = []
        for i in range(len(self.x)-2):
            x1 = self.x[i]
            x2 = self.x[i+1]
            x3 = self.x[i+2]
            z1 = self.z[i]
            z2 = self.z[i+1]
            z3 = self.z[i+2]
            a = np.sqrt((x1-x2)**2 + (z1-z2)**2)
            b = np.sqrt((x2-x3)**2 + (z2-z3)**2)
            c = np.sqrt((x3-x1)**2 + (z3-z1)**2)
            # Compute inverse radius of circle using surface of triangle (for which Heron's formula is used)
            k = np.sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))) / 4    # Heron's formula for triangle's surface
            den = a*b*c  # Denumerator; make sure there is no division by zero.
            if den == 0.0:  # Very unlikely, but just to be sure
                curvature.append(0.0)
            else:
                curvature.append(4*k / den)
        
        return np.array(curvature)
    
    def plotCurvature(self):
        curvature = self.getCurvature()
        plt.plot(np.linspace(0,0.6,100), -curvature)
        plt.xlabel('Position along the track')
        plt.ylabel('Curvature')
        plt.title('Curvature along the track')
        plt.grid(True)
        plt.show()
    
    def exportCurvature(self, filename):
        curvature = self.getCurvature()
        pd.DataFrame(curvature).to_csv(filename, index=False, header=False)
        print('Curvature exported to ' + filename)