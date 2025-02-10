import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np
materials = ['SUS304', 'NiTi_fixed', 'SiliconeRubber']

susList = []
nitiList = []
siliList = []

def getInfoFromPath(path):
    pathname = path.split('\\')[-1].split(".csv")[0].split('_')
    name = pathname[0]
    fratio = int(pathname[-2].split('fratio')[-1])
    return name, fratio

def average(list):
    sum = 0
    for val in list:
        sum += val
    return sum / len(list)

def y_intercept(x1, x2, y1, y2):
    return y1 - x1*(y2-y1)/(x2-x1)


folderpath = f'Python/curvaturedatas/*.csv'
fratios = [5, 7, 9, 11, 13, 15]
file_list = glob.glob(folderpath)
for force in fratios:
    for path in file_list:
        name, fratio = getInfoFromPath(path)
        print(f'reading {name}, {fratio}')
        if fratio == force:
            df = pd.read_csv(path)
            plt.plot(df['LRratio'], df['MaxCurv'],
                     marker='o' if name == 'SUS304' else 'x',
                     markersize = 10 if name == 'NiTi' else 5,
                     alpha = 0.5,
                     label = name + " Max.")
            plt.plot(df['LRratio'], df['MinCurv'],
                     marker='o' if name == 'SUS304' else 'x',
                     markersize = 10 if name == 'NiTi' else 5,
                     alpha = 0.5,
                     label = name + " Min.")
            rationd = np.array(df['LRratio'])[2:20]
            maxnd = np.array(df['MaxCurv'])[2:20]
            minnd = np.array(df['MinCurv'])[2:20]
            maxpoly = np.poly1d(np.polyfit(rationd, maxnd, 1))
            minpoly = np.poly1d(np.polyfit(rationd, minnd, 1))
            trendx = np.linspace(0,50,5)

            plt.plot(trendx, maxpoly(trendx), label = name + " Max. trendline")
            plt.plot(trendx, minpoly(trendx), label = name + " Min. trendline")

            y_int = y_intercept(df['LRratio'][2],df['LRratio'][15],df['MaxCurv'][2], df['MaxCurv'][15])
            avg = average(df['MinCurv'])
            print(f'{name}, y_intercept {y_int}, Lwr. avg {avg}')
    plt.legend()
    plt.title(f'Initial and Final curvature, fratio {force}')
    plt.xlabel('L/r')
    plt.ylabel('Curvature')
    plt.xlim((0,50))
    plt.tight_layout()
    plt.savefig(f'Python/figs/Curvature_fratio{force}.png')
    plt.show()
