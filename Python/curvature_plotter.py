from CosseratBeamFileHandle import CosseratBeamFile as cbf
import glob
import matplotlib.pyplot as plt
import numpy as np

path = 'Matlab/code/csvfiles/NiTi/*.csv'
file_list = glob.glob(path)

cbf_list = []
for file in file_list:
    cbf_list.append(cbf(file))

cbf_list.sort(key=lambda x: x.length)
intercepts = []

plt.figure()
plt.subplot(231)
tmplist = []
for beam in cbf_list:
    if beam.length == '0.10':
        plt.plot(beam.x, beam.z)
        plt.axis('equal')
        tmplist.append(float(beam.returnCurvLinIntercept()))
intercepts.append(tmplist)
plt.subplot(232)
tmplist = []
for beam in cbf_list:
    if beam.length == '0.20':
        plt.plot(beam.x, beam.z)
        plt.axis('equal')
        tmplist.append(float(beam.returnCurvLinIntercept()))
intercepts.append(tmplist)
plt.subplot(233)
tmplist = []
for beam in cbf_list:
    if beam.length == '0.30':
        plt.plot(beam.x, beam.z)
        plt.axis('equal')
        tmplist.append(float(beam.returnCurvLinIntercept()))
intercepts.append(tmplist)
plt.subplot(234)
for beam in cbf_list:    
    if beam.length == '0.10':
        x, z = beam.returnCurvature()
        plt.plot(x, z)
plt.subplot(235)
for beam in cbf_list:
    if beam.length == '0.20':
        x, z = beam.returnCurvature()
        plt.plot(x, z)
plt.subplot(236)
for beam in cbf_list:
    if beam.length == '0.30':
        x, z = beam.returnCurvature()
        plt.plot(x, z)
plt.show()
    
for row in intercepts:
    print(row)