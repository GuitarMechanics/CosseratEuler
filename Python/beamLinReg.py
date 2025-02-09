from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os
import cmd
import pandas as pd

# os.system('cls' if os.name == 'nt' else 'clear')
# os.system('cls')
materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']
# materialIndex = int(input('Select material: SUS304[0], NiTi[1], SiliconeRubber[2]'))
materialIndex = 1
# path = f'Matlab/code/csvfiles/{materials[materialIndex]}/*.csv'

for mat in materials:
    path = f'Matlab/code/csvfiles/{mat}/*.csv'
    file_list = glob.glob(path)
    cbf_list = []
    LRratioList = []
    curvFinalList = []
    curvInitiList = []
    curvCCEquList = []
    forceratio = [5, 7, 9, 11, 13, 15]
    for file in file_list:
        cbf_list.append(cbf(file))
        if cbf_list[-1].LRratio not in LRratioList:
            LRratioList.append(cbf_list[-1].LRratio)
        cbf_list.sort(key = lambda x: x.LRratio)
    for force in forceratio:
        length = []
        radius = []
        lr_ratios = []
        force_ratios = []
        max_curvatures = []
        eqcc_curvatures = []
        min_curvatures = []
        for beam in cbf_list:
            if beam.forceratio == force:
                # plotx = []
                # plotz = []
                # for xval in beam.x:
                #     plotx.append(xval / beam.length)
                # for zval in beam.z:
                #     plotz.append(zval / beam.length)
                # plt.plot(plotx, plotz)
                curvinfo = beam.getCurvatureInfo(normalized = True)
                lr_ratios.append(beam.LRratio)
                max_curvatures.append(curvinfo[0])
                min_curvatures.append(curvinfo[1])
                eqcc_curvatures.append(curvinfo[2])
                length.append(beam.length)
                radius.append(beam.radius)

                # plt.plot(beam.LRratio, curvinfo[0], marker='o', color = 'blue', label = 'max')
                # plt.plot(beam.LRratio, curvinfo[0], marker = 'o', color = 'red')
                # plt.plot(beam.LRratio, curvinfo[1], marker = 'o', color = 'green')
            dfdata = {'length' : length,
                      'radius' : radius,
                      'LRratio': lr_ratios,
                      'MaxCurv': max_curvatures,
                      'MinCurv': min_curvatures,
                      'EqCCval': eqcc_curvatures}
            pd.DataFrame(dfdata).to_csv(f'Python/{mat}_fratio{force}_curvaturedata.csv',index=None)