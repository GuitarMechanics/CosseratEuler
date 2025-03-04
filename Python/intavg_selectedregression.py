import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import clipboard
from math import pi
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
clipboard.copy
df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'curvature_reginfos', engine='openpyxl')
fratio = [5, 7, 9, 11, 13, 15]
for frat in fratio:
    LRratio = []
    interavg = []
    for i, row in df.iterrows():
        if row['Fratio'] == frat:
            LRratio.append(row['lrratio'])
            interavg.append(row['KR/AR'])
    LRfit = np.array(LRratio).reshape(-1,1)
    interavg = np.array(interavg).reshape(-1,1)
    model = LinearRegression().fit(LRfit, interavg)
    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = r2_score(interavg, model.predict(LRfit))
    print(f'Fratio: {frat} | Slope: {slope} | Intercept: {intercept} | R2 Score: {r2} | slope/FRatio: {slope/frat}')