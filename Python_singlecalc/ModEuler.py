import numpy as np

class ModEuler():
    def __init__(self, l, r):
        self.len = l
        self.rad = r
    
    def getInitCurv(self, tdl):
        return (tdl * (2*self.len + tdl) / (2 * (self.len * self.rad * (self.len - tdl))))
    
    def getTipAng(self, tdl):
        return tdl / self.rad
    
    def getTipPos(self, tdl, integral_resolution = 1000):
        segs = np.linspace(0,self.len, integral_resolution)
        k1 = self.getInitCurv(tdl)
        ka = tdl / (self.len * self.rad)
        ang_vals = []
        for val in segs:
            cval = (ka - k1) * val**2 / self.len + k1 * val
            ang_vals.append(cval)

        cvndarr = np.array(ang_vals)
        
        horpos = np.trapz(np.sin(cvndarr), segs)
        verpos = np.trapz(np.cos(cvndarr), segs)

        return [horpos, verpos]
