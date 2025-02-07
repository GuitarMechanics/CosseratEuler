from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os
import cmd

materials = ['SUS304', 'NiTi_past', 'SiliconeRubber']
materialIndex = int(input('Select material: SUS304[0], NiTi[1], SiliconeRubber[2]'))
path = f'Matlab/code/csvfiles/{materials[materialIndex]}/*.csv'

file_list = glob.glob(path)
cbf_list = []
LRratioList = []
curvFinalList = []
curvInitiList = []
curvCCEquList = []
for file in file_list:
    cbf_list.append(cbf(file))
    if cbf_list[-1].LRratio not in LRratioList:
        LRratioList.append(cbf_list[-1].LRratio)
    cbf_list.sort(key = lambda x: x.LRratio)
for beam in cbf_list:
    plt.plot(beam.LRratio, beam.getCurvatureInfo()[3])
plt.show()