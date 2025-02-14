import glob
from sklearn.metrics import r2_score
import pandas as pd
from CosseratBeamFileHandle import CosseratBeamFile as cbf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']

materialList = []
lengthList = []
rList = []
r2scoreList = []
forcelist = []

for material in materials:
    path = f'Matlab/code/csvfiles/{material}/*.csv'
    fileList = glob.glob(path)
    cbf_list = []
    for file in fileList:
        beam = cbf(file)
        materialList.append(beam.material)
        lengthList.append(beam.length)
        rList.append(beam.radius)
        forcelist.append(beam.forceratio)
        seg, curv = beam.getCurvature(normalize=False)
        fit_line = np.polyfit(seg, curv, 1)
        est_curv = seg * fit_line[0] + fit_line[1]
        r2 = r2_score(curv, est_curv)
        r2scoreList.append(r2)
        print(f'{beam.material}, {beam.length}, {beam.radius}, {r2}')

dfdata = {'Material' : materialList,
          'Length': lengthList,
          'radius': rList,
          'Fratio': forcelist,
          'R2score': r2scoreList}
df = pd.DataFrame(dfdata)

# pd.DataFrame(dfdata).to_csv('Python/curvature_r2.csv',index=None)