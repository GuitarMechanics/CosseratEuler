import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'autodetector_simple', engine='openpyxl')

r2_threshold = 0.70
related_columns = []

# 각 열 쌍에 대해 선형 회귀 수행 및 결정계수 계산
num_cols = df.shape[1]
columns = df.columns
m_col = df['intercept']
m_val = m_col.values

for i in range(num_cols):
    for j in range(i + 1, num_cols):  # 자기 자신 제외
        print(f"Processing columns {i} and {j}...")
        x1 = df[columns[i]]
        x2 = df[columns[j]]

        plt.figure(figsize=(15,10))
        plt.suptitle(f"{columns[i]}, {columns[j]} - m")
        plt.subplot(2,3,1)
        sns.scatterplot(x = x1, 
                        y = m_col, 
                        hue = df['Fratio'], 
                        palette = 'Set2')
        plt.title(f'{columns[i]} vs inter')

        plt.subplot(2,3,2)
        sns.scatterplot(x = x2,
                        y = m_col,
                        hue = df['Fratio'],
                        palette='Set2')
        plt.title(f'{columns[j]} vs inter')

        plt.subplot(2,3,3)
        xmul = x1*x2
        sns.scatterplot(x = xmul,
                        y = m_col,
                        hue = df['Fratio'],
                        palette = 'Set2')
        plt.title(f'{columns[i]}*{columns[j]} vs inter')

        plt.subplot(2,3,4)
        xdiv = x1/x2
        sns.scatterplot(x = xdiv,
                        y = m_col,
                        hue = df['Fratio'],
                        palette = 'Set2')
        plt.title(f'{columns[i]}/{columns[j]} vs inter')

        plt.subplot(2,3,5)
        xdiv2 = x2/x1
        sns.scatterplot(x = xdiv2,
                        y = m_col,
                        hue = df['Fratio'],
                        palette = 'Set2')
        plt.title(f'{columns[j]}/{columns[i]} vs inter')

        plt.savefig(f"Python/cor_detectorfig/inter_simpledataset_{i}_{j}.png")

# 결과 데이터프레임 생성
