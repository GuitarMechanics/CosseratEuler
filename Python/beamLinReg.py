from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os
import cmd

# os.system('cls' if os.name == 'nt' else 'clear')
# os.system('cls')
materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']
# materialIndex = int(input('Select material: SUS304[0], NiTi[1], SiliconeRubber[2]'))
materialIndex = 1
path = f'Matlab/code/csvfiles/{materials[materialIndex]}/*.csv'

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
            plt.plot(beam.LRratio, curvinfo[0], marker='o', color = 'blue', label = 'max')
            # plt.plot(beam.LRratio, curvinfo[0], marker = 'o', color = 'red')
            # plt.plot(beam.LRratio, curvinfo[1], marker = 'o', color = 'green')
    plt.title(f'Force ratio = {force}')
    plt.show()