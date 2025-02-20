import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_excel('Python/curvature_reginfo_xlvers.xlsx',sheet_name = 'curvature_reginfos', engine='openpyxl')

sns.scatterplot(
    data=df,
    x = 'lrratio',
    y = 'inter/avg',
    hue = 'Fratio',
    palette = 'Set2',
)
plt.show()