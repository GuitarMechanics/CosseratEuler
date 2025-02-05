from CosseratBeamFileHandle import CosseratBeamFile as cbf
import matplotlib.pyplot as plt
import glob
import os

os.system('cls')

path = 'Matlab/code/csvfiles/SiliconeRubber/*.csv'
file_list = glob.glob(path)
cbf_list = []
for file in file_list:
    cbf_list.append(cbf(file))
# cbf_list.sort(key=lambda x: x.length)
cbf_list.sort(key = lambda x: x.LRratio)

dimensions = []
for beam in cbf_list:
    if (beam.radius, beam.length, beam.LRratio) not in dimensions:
        dimensions.append((beam.radius, beam.length, beam.LRratio))
dimensions.sort(key = lambda x:x[2])
for dim in dimensions:
    ratiolist = []
    maxcurvaturelist = []
    mincurvaturelist = []
    difflist = []
    for beam in cbf_list:
        if beam.forceratio == 15:
            ratiolist.append(beam.LRratio)
            _, curvature = beam.returnCurvature(normalized = True)
            maxcurvaturelist.append(curvature[0])
            mincurvaturelist.append(curvature[-1])
            difflist.append(curvature[0] - curvature[-1])
plt.figure()
plt.plot(ratiolist, maxcurvaturelist, marker='o', label = "max")
plt.plot(ratiolist, mincurvaturelist, marker='o', label = "min")
plt.plot(ratiolist, difflist, marker='o', label = "diff")
plt.legend()
plt.tight_layout()
plt.show()


exit()
for dim in dimensions:
    tmplist = []
    for beam in cbf_list:
        if (beam.radius, beam.length) == (dim[0], dim[1]):
            tmplist.append(beam)
    tmplist.sort(key = lambda x:x.forceratio)
    plt.figure()
    for beam in tmplist:
        plt.subplot(131)
        plt.plot(beam.x, beam.z, label = str(beam.forceratio))
        plt.title("Configuration")
        plt.subplot(132)
        segment, curv = beam.returnCurvature(normalized = True)
        plt.plot(segment, curv, label=str(beam.forceratio))
        plt.title("Curvature")
        plt.subplot(133)
        segment, ang = beam.getAngles(normalized = False)
        plt.plot(segment, ang, label = str(beam.forceratio))
    plt.legend()
    # plt.axis('equal')
    plt.title(f'Radius {dim[0]}, Length {dim[1]}, LR_ratio = {dim[2]}')
    # plt.ylim(0,0.15)
    plt.tight_layout()
    plt.show()
    # dummy = input('Draw next?')