import ModEuler as me
import numpy as np
import pandas as pd

std_tdls = np.linspace(-8.75, 8.75, 15)
off_tdls = np.linspace(-8.84, 8.84, 21)

tdcr = me.ModEuler(61,3.4)

std_posx = []
std_posy = []
off_posx = []
off_posy = []

for tdl in std_tdls:
    hpos, vpos = tdcr.get_tippos(tdl)
    std_posx.append(hpos)
    std_posy.append(vpos)
for tdl in off_tdls:
    hpos, vpos = tdcr.get_tippos(tdl)
    off_posx.append(hpos)
    off_posy.append(vpos)

dfdata_std = {'std_tdl' : std_tdls,
          'std_x' : std_posx,
          'std_y' : std_posy}

dfdata_off = {'off_tdl' : off_tdls,
              'off_x' : off_posx,
              'off_y' : off_posy}

pd.DataFrame(dfdata_std).to_csv('Python/std_positions.csv',index=None)
pd.DataFrame(dfdata_off).to_csv('Python/off_positions.csv',index=None)