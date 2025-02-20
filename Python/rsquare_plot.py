import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('Python/curvature_r2_normalized.csv')
print(df)
plt.figure(figsize = (16,9))
plt.suptitle('R^2 score distribution')
plt.subplot(221)
plt.title('Overall')
sns.boxplot(y = 'Material',
            x = 'R2score',
            fill = False,
            data = df,
            palette='Set2').set(
                ylabel = "Material",
                xlabel = "R^2"
            )
plt.subplot(222)
plt.title('By the Length')
sns.boxplot(y = 'R2score',
            x = 'Length',
            fill = False,
            data = df,
            palette='Set2',
            hue = 'Material').set(
                xlabel = "length",
                ylabel = "R^2"
            )
plt.subplot(223)
plt.title('By the radius')
sns.boxplot(y = 'R2score',
            x = 'radius',
            fill = False,
            data = df,
            palette='Set2',
            hue = 'Material').set(
                xlabel = "radius",
                ylabel = "R^2"
            )
plt.subplot(224)
plt.title('By the force ratio')
sns.boxplot(y = 'R2score',
            x = 'Fratio',
            fill = False,
            data = df,
            palette='Set2',
            hue = 'Material').set(
                xlabel = "Force ratio",
                ylabel = "R^2"
            )

plt.tight_layout()
plt.savefig('Python/Curvature_distribution_normalized.png')
plt.show()
