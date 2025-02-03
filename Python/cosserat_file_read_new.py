from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import glob
path = 'Matlab/code/csvfiles/SUS304/*.csv'
file_list = glob.glob(path)
cbf_list = []
for file in file_list:
    cbf_list.append(cbf(file))
cbf_list.sort(key=lambda x: x.length)
plt.figure()
for beam in cbf_list:
    if float(beam.length) == 0.059:
        plt.plot(beam.x, beam.z, label = "r = " + beam.radius)
plt.axis('equal')
plt.legend()
plt.show()
