import CosseratBeam as cb
import numpy as np
import matplotlib.pyplot as plt
import FresnelEuler as fe
import pandas as pd
import statsmodels.api as sm
resolution = 50
length = 0.6
forceArr = [70,80,90,100,110,120]
def getFileName(force,resolution,length):
    return f"Matlab\code\csvfiles\SingleSectionCR_L{length:.2f}_N{resolution:0.0f}_r0.0100_Tt1{force:.2f}_Tt20.00.csv"

eulerCurve = fe.FresnelEuler(length,3)

beams = []
tdls = []
for force in forceArr:
    beams.append(cb.CosseratBeam(force,resolution,length))

beams_initcurvature = []
posArry = np.linspace(0,length,resolution)
plt.figure()
plt.subplot(1,2,1)
curvature_xinter = []
for beam in beams:
    plt.plot(np.linspace(0,length,resolution)[5:-5], beam.getCurvature()[5:-5])
    # pd.DataFrame(data={'Position': np.linspace(0,length,resolution),
    #                    'Curvature': beam.getCurvature()}).to_csv(f'curvature_{beam.force}.csv')
    X = sm.add_constant(np.linspace(0,length,resolution)[5:-5])
    model = sm.OLS(beam.getCurvature()[5:-5],X)
    results = model.fit()
    print(-results.params[0]/results.params[1])

pd.DataFrame(curvature_xinter).to_csv('curvature_xinter.csv',index=False,header=False)
newcoss = cb.CosseratBeam(200,100,0.6)
plt.plot(np.linspace(0,length,100)[5:-5],newcoss.getCurvature()[5:-5])

    # print(np.max(beam.getCurvature()))
    # beams_initcurvature.append(np.max(beam.getCurvature()[0]))
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

exit()
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