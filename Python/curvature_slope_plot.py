import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
df = pd.read_csv('Python/curvature_reginfos_centcurvfix.csv')

sns.regplot(data = df, x = 'AvgCurv', y = 'intercept',
            scatter_kws={"color" : "gray", "alpha" : 0.7},
            line_kws= {"color" : "black"},
            ci = 99)
plt.xlabel('Average curvature' + r"$(\kappa_{avg})$" + "  " + r"$(m^{-1})$")
plt.ylabel('Initial curvature'+ r"$(\kappa_{1})$"+ "  " + r"$(m^{-1})$")
# plt.text(60,20,r"$\kappa_{1} = 1.085\kappa_{avg} + 0.429$" + '\n' + r"$R^{2} = 0.997$")

x = np.array(df['AvgCurv'])
y = np.array(df['intercept'])
X = x[:, np.newaxis]
a, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
print(a)
plotx = np.array([0,110])
ploty = np.array([0,110*a[0]])
plt.plot(plotx, ploty)
plt.tight_layout()
# plt.savefig('Python/curvature_initavg_relationship.png')
plt.show()
exit()
exit()
pfit = np.polyfit(df['AvgCurv'], df['intercept'],1)
print(pfit[0], pfit[1])
polydraw = np.array([0,110])

plt.plot(df['AvgCurv'], df['intercept'], linewidth = 0, marker = 'o')
pfitxy = polydraw * pfit[0] + pfit[1]
plt.plot(polydraw, pfitxy)
plt.show()