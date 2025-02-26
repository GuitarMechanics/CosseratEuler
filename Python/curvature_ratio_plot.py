import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'curvature_reginfos', engine='openpyxl')
# plt.xlim([0,0.12])
# plt.ylim([0,0.012])
sns.scatterplot(
    data=df,
    x = df['inter/avg/tang*LR'],
    y = df['inter/avg-1'],
    hue = 'Fratio',
    palette = 'Set2'
    # hue = 'tipang'
)#.set_xlabel('1/LR')
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
