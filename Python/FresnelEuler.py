import os
import numpy as np
import pandas as pd
import math
from math import pi
import scipy as sp
from scipy.optimize import fsolve
class FresnelEuler:
    def __init__(self, armlen, diskrad):
        self.armLen = armlen
        self.r_val = diskrad
        
    def fresnel_sin(self, val):
        retval = 0
        for n in range(20):
            retval += ((-1)**n * (pi/2)**(2*n+1) * val**(4*n+3)) / ((4*n+3)*math.factorial(2*n+1))
        return retval
    def fresnel_cos(self, val):
        retval = 0
        for n in range(20):
            retval += ((-1)**n * (pi/2)**(2*n) * val**(4*n+1)) / ((4*n+1)*math.factorial(2*n))
        return retval

    def getDistPos(self, tdl):
        #tdl = float(input("TDL: "))
        if(tdl == 0):
            return ([0,self.armLen])
        a_val = 2*tdl / self.r_val
        sqrt_val = math.sqrt(math.pi / abs(a_val))
        # print("a_val ",a_val, "sqrt_val ",sqrt_val)

        xpos = sqrt_val * (self.fresnel_cos(1/sqrt_val) * math.sin(a_val/2) - \
                        self.fresnel_sin(1/sqrt_val)*math.cos(a_val/2))
        ypos = sqrt_val * (self.fresnel_cos(1/sqrt_val) * math.cos(a_val/2) + \
                        self.fresnel_sin(1/sqrt_val)*math.sin(a_val/2))
        xpos *= self.armLen
        ypos *= self.armLen

        # print(xpos,ypos)
        return ([xpos,ypos])

    def curvCalculator(self,tdl,pos):
        if(tdl == 0):
            return ([0,pos])
        norm_xval = pos/self.armLen
        a_val = 2*abs(tdl) / self.r_val
        sqrt_val = math.sqrt(math.pi / a_val)
        # print("a_val ",a_val, "sqrt_val ",sqrt_val)

        ypos = sqrt_val * (math.cos(a_val/2)*(self.fresnel_cos((norm_xval - 1)/sqrt_val) + self.fresnel_cos(1/sqrt_val)) + \
                        math.sin(a_val/2)*(self.fresnel_sin((norm_xval - 1)/sqrt_val) + self.fresnel_sin(1/sqrt_val)))
        xpos = sqrt_val * (math.sin(a_val/2)*(self.fresnel_cos((norm_xval - 1)/sqrt_val) + self.fresnel_cos(1/sqrt_val)) - \
                        math.cos(a_val/2)*(self.fresnel_sin((norm_xval - 1)/sqrt_val) + self.fresnel_sin(1/sqrt_val)))
        xpos *= self.armLen
        ypos *= self.armLen
        return([xpos,ypos])
    
    def getTDLfromInitCurvature(self,intialCurvature):
        return intialCurvature*self.r_val*self.armLen/2
    
    def curveCalculatorfromInitCurvature(self,intialCurvature, pos):
        tdl = self.getTDLfromInitCurvature(intialCurvature)
        return self.curvCalculator(tdl,pos)
    
    def getBestMatchTDL(self, tipx, tipz):
        tdlstart = 0.01
        tdl = fsolve(lambda x: np.hypot(\
            self.getDistPos(x)[0] - tipx,self.getDistPos(x)[1]-tipz), tdlstart)
        print(tdl)
        return tdl