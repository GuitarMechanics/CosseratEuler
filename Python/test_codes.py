# from CosseratBeamFileHandle import CosseratBeamFile as cbf
import glob
import numpy as np
# from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

r = np.linspace(1,10,19)
l = np.linspace(50,500,37)
ratios = []
for rval in r:
    for lval in l:
        if lval/rval >= 5 and lval/rval <= 100:
            ratios.append(lval/rval)
ratios.sort()
print(len(ratios))
plt.figure()
plt.plot(np.linspace(0,len(ratios),len(ratios)),ratios,marker='o')
plt.show()

exit()
# path = 'Matlab/code/csvfiles/NiTi/*.csv'
# file_list = glob.glob(path)
# cbf_list = []
# for file in file_list:
#     cbf_list.append(cbf(file))

# cbf_list.sort(key=lambda x: x.length)

# beam = cbf_list[0]
# pos, curvature = beam.returnCurvature()

# linearreg = LinearRegression().fit(np.array(curvature).reshape(-1,1),np.array(pos).reshape(-1,1))

# print(linearreg.intercept_[0])
# print(linearreg.coef_[0][0])


# linearreg2 = LinearRegression().fit(np.linspace(0,10,20).reshape(-1,1),(np.linspace(0,50,20)+1).reshape(-1,1))

# print(linearreg2.intercept_[0])