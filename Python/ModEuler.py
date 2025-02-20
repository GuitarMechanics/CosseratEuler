import os
import numpy as np
import pandas as pd
import math
from math import pi
import scipy as sp
from scipy.integrate import quad
class ModEuler:
    def __init__(self, armlen, diskrad):
        self.armLen = armlen
        self.r_val = diskrad
        
    def get_TDLfromAngle(self, angle):
        return self.r_val * angle

    def get_avgCurvature(self, tdl):
        return tdl / (self.armLen * self.r_val)
    
    def get_tippos(self, tdl):
        avgcurvature = self.get_avgCurvature(tdl)
        hor = lambda x: np.sin(-0.1*avgcurvature*x**2/self.armLen \
                               + 1.1*avgcurvature * x)
        ver = lambda x: np.cos(-0.1*avgcurvature*x**2/self.armLen \
                               + 1.1*avgcurvature * x)
        horpos = quad(hor,0,self.armLen)[0]
        verpos = quad(ver,0,self.armLen)[0]
        return horpos, verpos