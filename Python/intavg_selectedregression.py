import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import clipboard
from math import pi
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np

df = pd.read_excel('Python/curvature_reginfo_xlvers_gpt.xlsx',sheet_name = 'pade_sols', engine='openpyxl')
fratio = [5, 7, 9, 11, 13, 15]
for frat in fratio:
    xvals = []
    yvals = []
    for i, row in df.iterrows():
        if row['Fratio'] == frat:
            xvals.append(row['tang'])
            yvals.append(row['Lint+tang'])
    x = np.array(xvals).reshape(-1,1)
    y = np.array(yvals).reshape(-1,1)
    model = LinearRegression().fit(x, y)
    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = r2_score(y, model.predict(x))
    print(f'Fratio: {frat} | Slope: {slope} | Intercept: {intercept} | R2 Score: {r2} | slope/FRatio: {slope/frat} | X-inter: {-intercept/slope}')