import glob
import pandas as pd
from sklearn.metrics import r2_score
from CosseratBeamFileHandle import CosseratBeamFile as cbf
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']

materialList = []
lengthList = []
rList = []
forcelist = []
slopeList = []
interList = []
centrList = []
tdlList = []
curvDiffList = []

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
        slopeList.append(fit_line[0])
        interList.append(fit_line[1])
        centcurv = np.mean(curv)
        tdl = centcurv * beam.length * beam.radius
        centrList.append(centcurv)
        tdlList.append(tdl)
        curvDiffList.append(fit_line[1] - centcurv)
        print(f'Processing : {beam.material}, {beam.length}, {beam.radius}')

dfdata = {'Material' : materialList,
          'Length': lengthList,
          'radius': rList,
          'Fratio': forcelist,
          'TDL': tdlList,
          'slope': slopeList,
          'intercept': interList,
          'AvgCurv': centrList,
          'CurvDiff': curvDiffList
          }
df = pd.DataFrame(dfdata)

pd.DataFrame(dfdata).to_csv('Python/curvature_reginfos_centcurvfix.csv',index=None)