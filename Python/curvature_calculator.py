import CosseratBeam as cb
import numpy as np
import matplotlib.pyplot as plt
import FresnelEuler as fe
resolution = 50
length = 0.6
forceArr = [10,20,30,40,50,60]
def getFileName(force,resolution,length):
    return f"Matlab\code\csvfiles\SingleSectionCR_L{length:.2f}_N{resolution:0.0f}_r0.0100_Tt1{force:.2f}_Tt20.00.csv"

eulerCurve = fe.FresnelEuler(length,3)

beams = []
tdls = []
for force in forceArr:
    beams.append(cb.CosseratBeam(getFileName(force,resolution,length)))

beams_initcurvature = []
posArry = np.linspace(0,length,resolution)
plt.figure()
plt.subplot(1,2,1)
for beam in beams:
    plt.plot(np.linspace(0,length,resolution), beam.getCurvature())
    print(np.max(beam.getCurvature()))
    beams_initcurvature.append(np.max(beam.getCurvature()[0]))
plt.xlabel('Position along the track')
plt.ylabel('Curvature')
plt.title('Curvature along the track')

plt.subplot(1,2,2)
plt.axis('equal')
for beam in beams:
    plt.plot(beam.x, beam.z)
    tdl = eulerCurve.getBestMatchTDL(beam.x[-1],beam.z[-1])
    tdls.append(tdl)
    eulerx = []
    eulerz = []
    for i in range(resolution):
        eulerx.append(eulerCurve.curvCalculator(tdl,posArry[i])[0])
        eulerz.append(eulerCurve.curvCalculator(tdl,posArry[i])[1])
    plt.plot(eulerx,eulerz,linestyle='--')
plt.xlabel('x')
plt.ylabel('z')
plt.title('Configurations')
plt.show()

plt.figure()
beamindex = 0
for beam in beams:
    poserr = []
    for i in range(resolution):
        eulerpos = eulerCurve.curvCalculator(tdls[beamindex],posArry[i])
        poserr.append(np.hypot(beam.x[i]-eulerpos[0],beam.z[i]-eulerpos[1]) / length * 100)
    plt.plot(posArry,poserr)
    beamindex += 1
plt.show()