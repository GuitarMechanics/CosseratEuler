from ModEuler import ModEuler as me
import os

os.system('cls') if os.name == 'nt' else os.system('clear')
tdlarrs = [-8.75]

while tdlarrs[-1] < 8.75:
    tdlarrs.append(tdlarrs[-1] + 1.25)


tdcr = me(62, 3.4)

for val in tdlarrs:
    print(tdcr.getTipPos(tdl = val))