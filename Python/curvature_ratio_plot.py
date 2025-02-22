import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'curvature_reginfos', engine='openpyxl')

sns.scatterplot(
    data=df,
    x = 'inter/avg',
    y = 'm',
    # hue = 'Fratio',
    # palette = 'Set2'
    hue = 'lrratio'
)
plt.show()

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