import math
import numpy as np

print(np.linspace(0,0.1,50))

exit()

filepath = "Matlab/code/csvfiles/NiTi/NiTi_L0.20_N50_r0.0100_Tt143982.30.csv"
filename = filepath.split("/")[-1].split(".csv")[0]
print(filename)

filename = filename.split("_")
forceval = filename[4].split("Tt1")[1]
length = filename[1].split("L")[1]
resolution = filename[2].split("N")[1]
print(forceval)
print(length)
print(resolution)
exit()
filename = "NiTi_L0.10_N50_r0.0100_Tt1285884.93"



print(filename)
if filename[0] == "NiTi":
    print("NiTi")
    E_val = 70e9
    A_val = 0.01**2*math.pi
    forceval = filename[4].split("Tt1")[1]
    print(np.round(float(forceval)/(E_val*A_val)*1e3))
    forceval = int(np.round(float(forceval)/(E_val*A_val)*1e3))*1e-3
    print(forceval)