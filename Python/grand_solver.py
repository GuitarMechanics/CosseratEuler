import sympy as sp
import os

os.system('cls') if os.name == 'nt' else os.system('clear')

# Define the variables

x, k1, r, Lb, dt = sp.symbols('x \\kappa_1 r L t')


ps = dt * (2 * Lb + dt) / (2 * Lb * r * (Lb - dt))
ka = dt/(Lb*r)

# ps, ka = sp.symbols("\\kappa_p \\kappa_a")

tsteqn = 2*(ps/ka+k1/ps) + ((ps + k1) / (ps + ka) + (ps + ka) / (ps + k1))/2 -2 - k1/ka
sp.pprint(sp.simplify(tsteqn))
sol = sp.solve(tsteqn,k1)

for val in sol:
    print('sol')
    sp.pprint(sp.latex(sp.simplify(val)))


exit()

curvat = (ka - k1)*x / (Lb/2) + k1

print(sp.latex(sp.simplify(curvat)))

ang = sp.integrate(curvat, x)
print(sp.latex(sp.simplify(ang)))

exit()
eqn = (ps / ka - k1 / ka)/(k1 / ps - 1) + ka**2/ps
print('EQUATION')
print(sp.latex(sp.simplify(eqn)))
print('------------------------')

sol = sp.solve(eqn, k1)
print('SOLUTION k1')
for val in sol:
    print(sp.latex(sp.simplify(val)))
exit()

u = 2 * r * (k1 - ka)

# eqn = (2*(Lb - dt) * r * (k1 - ka)) / (Lb * (1 - k1 * r)) + 1 
eqn = 1 + (2 * r * (Lb - dt) * (k1 - ka)) / (Lb*(1 - k1 * r))

print(sp.latex(sp.simplify(eqn)))

pade_11 = (2 + u) / (2 - u)

pade_21 = (1 + 2 * u / 3 + u**2 / 6) / (1 - u / 3)
pade_12 = (1 + u / 3) / (1 - 2 * u / 3 + u**2 / 6)
pade_22 = (1 + u / 2 + u**2 / 12) / (1 - u / 2 + u**2 / 12)

pade_13 = (u / 4 + 1) / (-u**3 / 24 + u**2 / 4 -3*u / 4 + 1)
pade_23 = (u**2 / 20 + 2* u / 5 + 1) / (-u**3 / 60 + 3*u**2 / 20 - 3*u / 5 + 1)
pade_31 = (u**3 / 24 + u**2 / 4 + 3*u / 4 + 1) / (-u / 4 + 1)
pade_32 = (u**3 / 60 + 3*u**2 / 20 + 3*u / 5 + 1) / (u**2 / 20 - 2*u / 5 + 1)
pade_33 = (1 + u / 2 + u**2 / 10 + u**3 / 120) / (1 - u / 2 + u**2 / 10 - u**3 / 120)



pade_list = [pade_11, pade_21, pade_12, pade_22, pade_13, pade_23, pade_31, pade_32, pade_33]
pade_name = ['pade_11', 'pade_21', 'pade_12', 'pade_22', 'pade_13', 'pade_23', 'pade_31', 'pade_32', 'pade_33']
nameindex = 0
for pade in pade_list:
    print("Solve for (u)" + pade_name[nameindex])
    nameindex += 1
    print("-------------------------------------------------")
    sol = sp.solve(eqn - pade, ka)
    for val in sol:
        print(sp.latex(sp.simplify(val)))
        print("-------------------------------------------------")
    print("=================================================")
exit()

kt_subs = sp.Subs(kt, (ka), (dt/(Lb*r)))
print("kt------")
print(sp.latex(sp.simplify(kt)))
# print(sp.latex(sp.simplify(kt_subs)))
# Perform the definite integration with symbolic limits
integral = sp.simplify(sp.integrate(sp.simplify(kt), (x, 0, Lb-dt)))

# Output the result in LaTeX format
int_equation = integral - dt/r

latex_output = sp.latex(sp.simplify(integral))

print("integral------")
print(latex_output)


print("int_equation = 0 --- == -")
print(sp.latex(sp.simplify(int_equation)))

print("int_equation, subs kavg as t/Lr")
print(sp.latex(sp.simplify(sp.Subs(int_equation, (ka), (dt/(Lb*r))))))

eqn_simplever = sp.exp(2*r*(k1-ka)) - 1 - \
    (2*r*(Lb-dt)*(k1-ka) / Lb*(1-r*k1))

print("eqn_simplever = 0 --- == -")
print(sp.latex(eqn_simplever))

eqn_solve = sp.solve(eqn_simplever, k1, dict = True)

print("eqn_solve")
sp.pprint(eqn_solve)
exit()