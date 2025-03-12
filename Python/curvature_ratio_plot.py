import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import clipboard
from math import pi
import numpy as np
clipboard.copy('C:\\Users\\user\\Desktop\\codes\\cosseratEuler\\Python\\excelfigs')
df = pd.read_excel('Python/curvature_reginfo_xlvers_gpt.xlsx',sheet_name = 'pade_sols', engine='openpyxl')
# plt.xlim([0,52])
# plt.ylim([0,1.2])
sns.scatterplot(
    data=df,
    x = df['lrratio'],
    # y = (df['tipang'])/df['inter/avg'],
    y = df['mod2_pct'],
    hue = 'Fratio',
    palette = 'Set2'
    # hue = 'tipang'
)#.set_xlabel('LR/(inter/avg)')
# plt.title('')
plt.show()
#path to copy C:\Users\user\Desktop\codes\cosseratEuler\Python\excelfigs

# # Create a polar plot
# plt.figure()
# ax = plt.subplot(111, projection='polar')
# sns.scatterplot(
#     data=df,
#     x='slope(atan)deg',
#     y='lrratio',
#     ax=ax
# )
# plt.show()
