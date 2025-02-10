import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

class CosseratBeamFile():
    def __init__(self, filelocation):
        '''
        If the filelocation is given, the data is read from the file.
        Also reveals the beam's information from the file name.
        The option normalized peforms re-scaling the size of the beam. After normailzation, the beam's length becomes 1.
        '''
        filename = filelocation.split("\\")[-1].split(".csv")[0]
        filename = filename.split("_")

        if filename[0] == "NiTi":
            self.E_val = 70e9
        elif filename[0] == "SUS304":
            self.E_val = 200e9
        elif filename[0] == "SiliconeRubber":
            self.E_val = 50e6
        else:
            raise Exception("The material is not recognized.")
        self.material = filename[0]
        forceval = filename[4].split("Tt1")[1]
        self.forceval = float(forceval)
        self.length = float(filename[1].split("L")[1])
        self.resolution = int(filename[2].split("N")[1])
        self.radius = float(filename[3].split("r")[1])
        A_val = self.radius**2*np.pi
        I_val = self.radius**4*np.pi/64
        self.forceratio = int(np.round(self.forceval*self.length*self.radius/(self.E_val*I_val)))
        self.LRratio = np.round(self.length / self.radius, 3)

        self.segments = np.linspace(0,self.length,self.resolution)
        
        self.df = pd.read_csv(filelocation,header=None)
        self.x = []
        self.z = []
        for i, row in self.df.iterrows():
            self.x.append(row[0])
            self.z.append(row[2])

    def getAngles(self, normalized = False):
        retList = []
        for i in range(self.resolution-1):
            retList.append(90-np.rad2deg(np.arctan2(self.z[i+1]-self.z[i],\
                                      self.x[i+1]-self.x[i])))
        if normalized:
            return self.segments[:-1], np.array(retList)
        else:
            return self.segments[:-1], np.array(retList)
      
    def getCurvature(self,normalize = False):
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
        if not normalize:
            return self.segments[2:-2], curvature[2:-2]
        else:
            return self.segments[2:-2] / self.length, curvature[2:-2] * self.length


    def returnCurvature(self, normalized = False):
        segment, curvature = self.getCurvature()
        if normalized == False:
            return segment, curvature
        else:
            return np.linspace(0,1,self.resolution)[2:-2], curvature * self.radius**2 / self.length
    
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
    
    def getEquivCCcurvature(self, normalized = False):
        angle = np.deg2rad(self.getAngles()[-1][-1])
        if normalized:
            return angle
        else:
            return angle/self.length

    def getCurvatureInfo(self, normalized = False, exportBeamInfo = False):
        '''
        returns: [initial curvature, final curvature, equivalent CC curvature, lr_ratio, forceratio]
        If exportBeamInfo == True, length, radius are appended
        '''
        seg, curvature = self.getCurvature(normalized)
        X = seg[:10].reshape(-1,1)
        curvature = curvature[:10]

        model = LinearRegression()
        model.fit(X,curvature)
        curv_initfinal = model.predict(np.array([[0], [1 if normalized else self.length]]))
        curv_equivcc = self.getEquivCCcurvature(normalized)

        retList = [curv_initfinal[0], curv_initfinal[1], curv_equivcc, self.LRratio, self.forceratio]
        if exportBeamInfo==False:
            return retList
        else:
            retList.append(self.length)
            retList.append(self.radius)
            return retList
