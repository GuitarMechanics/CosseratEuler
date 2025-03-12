import sympy as sp
from sympy.solvers.solveset import solveset_real
from sympy.solvers.solvers import solve
from scipy.optimize import fsolve
# Define the variables

x, k1, ka, r, Lb, dt = sp.symbols('x \\kappa_1 \\kappa_avg r L t')

ar = (1 + k1/ka)/4
kr = k1/ka
eqn = 16*ar + kr - 8 - 1/(ar*kr)

# eqn = eqn.subs({ka:dt/(Lb*r)})
solveval = sp.solve(eqn, k1)
for val in solveval:
    print(sp.latex(val))
    print('------------------')

exit()
a, b = sp.symbols('a b')

fx = sp.exp(2*r*(k1-ka)) / (k1-ka) -1/(k1-ka) -2*r*(Lb-dt)*(1-k1*r)/Lb
fx = fx.subs({ka:dt/(Lb*r)})

solveeqn = sp.lambdify(k1, fx.subs({Lb: 0.01, r:0.001, dt:0.001071}),modules=['numpy'])
solveval = fsolve(solveeqn, 0.001071/(0.01*0.001) * 2, xtol = 1e-6)
print(solveval)
# print(sp.latex(sp.simplify(solve[0])))