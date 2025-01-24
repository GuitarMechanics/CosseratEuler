import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

class CosseratBeamFile():
    def __init__(self, filelocation):
        '''
        If the filelocation is given, the data is read from the file.
        Also reveals the beam's information from the file name.

        CAUTION: The forceratio is scaled by 1e-3.
        example: 123 -> 0.123
        '''
        filename = filelocation.split("/")[-1].split(".csv")[0]
        filename = filename.split("_")

        if filename[0] == "NiTi":
            E_val = 70e9
        elif filename[0] == "SpringSteel":
            E_val = 207e9
        else:
            raise Exception("The material is not recognized.")
        
        A_val = 0.01**2*np.pi
        forceval = filename[4].split("Tt1")[1]
        forceval = int(np.round(float(forceval)/(E_val*A_val)*1e3))
        self.forceratio = forceval
        self.length = filename[1].split("L")[1]
        self.resolution = filename[2].split("N")[1]
        
        self.df = pd.read_csv(filelocation,header=None)
        self.x = []
        self.z = []
        for i, row in self.df.iterrows():
            self.x.append(row[0])
            self.z.append(row[2])

    # def __init__(self, force,resolution,length):
    #     '''
    #     If the force, resolution and length are given, the data is read from the file. Default location is Matlab\code\csvfiles, starting with SingleSectionCR_L.
    #     '''
    #     self.df = pd.read_csv(f"Matlab\code\csvfiles\SingleSectionCR_L{length:.2f}_N{resolution:0.0f}_r0.0100_Tt1{force:.2f}_Tt20.00.csv",header=None)
    #     self.x = []
    #     self.z = []
    #     self.force = force
    #     self.resolution = resolution
    #     self.length = length
    #     for i, row in self.df.iterrows():
    #         self.x.append(row[0])
    #         self.z.append(row[2])
        
    def getCurvature(self,returnIntercept = False):
        """
        Get the curvature of the track.
        if returnIntercept is True, the linear regression intercept of the curvature is returned.
        """
        x = self.x
        z = self.z

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
        if not returnIntercept:
            return curvature
        else:
            pos = np.linspace(0,float(self.length),int(self.resolution))[-2:2]
            linearreg = LinearRegression().fit(np.array(curvature)[2:-2].reshape(-1,1),np.array(pos).reshape(-1,1))
            return curvature, linearreg.intercept_[0]


    def returnCurvature(self, returnIntercept = False):
        curvature = self.getCurvature()
        return np.linspace(0,float(self.length),int(self.resolution))[2:-2], curvature[2:-2]
    
    def returnCurvLinIntercept(self):
        '''
        Returns the linear regression intercept of the curvature.
        '''
        pos, curvature = self.returnCurvature()
        linearreg = LinearRegression().fit(np.array(curvature).reshape(-1,1),np.array(pos).reshape(-1,1))
        return linearreg.intercept_[0]
    
    def returnConfiguration(self):
        return self.x, self.z
    
    def plotCurvature(self):
        curvature = self.getCurvature()
        plt.plot(np.linspace(0,self.length,self.resolution), -curvature)
        plt.xlabel('Position along the track')
        plt.ylabel('Curvature')
        plt.title('Curvature along the track')
        plt.grid(True)
        plt.show()
    
    def exportCurvature(self, filename):
        curvature = self.getCurvature()
        pd.DataFrame(curvature).to_csv(filename, index=False, header=False)
        print('Curvature exported to ' + filename)